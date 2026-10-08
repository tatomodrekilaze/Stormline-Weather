#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr std::size_t kSites = 20;
constexpr double kPi = 3.14159265358979323846;

template <typename Number>
struct State {
    std::vector<Number> circulation;
    std::vector<Number> humidity;
    std::vector<Number> cloud;
    std::vector<Number> rain;
};

template <typename Number>
struct Rates {
    State<Number> change;
};

template <typename Number>
struct RunResult {
    State<Number> finalState;
    std::vector<Number> rainBySite;
    std::vector<bool> stormBySite;
};

struct Options {
    std::size_t steps = 6000;
    double timeStep = 0.01;
    std::string outputPath = "stormline.csv";
    std::string svgPath = "stormline.svg";
};

struct SensitivityAtlas {
    std::array<std::array<double, kSites>, kSites> stateDifference{};
    std::array<std::array<bool, kSites>, kSites> stormStatusChanged{};
    double largestStateDifference = 0.0;
};

std::size_t wrapIndex(long index) {
    const long count = static_cast<long>(kSites);
    return static_cast<std::size_t>((index % count + count) % count);
}

template <typename Number>
State<Number> makeInitialState() {
    State<Number> state{
        std::vector<Number>(kSites),
        std::vector<Number>(kSites),
        std::vector<Number>(kSites, Number{0}),
        std::vector<Number>(kSites, Number{0}),
    };

    const Number forcing = static_cast<Number>(8.0);
    const Number humidityBase = static_cast<Number>(0.82);
    const Number humidityRange = static_cast<Number>(0.12);

    for (std::size_t i = 0; i < kSites; ++i) {
        const Number position = static_cast<Number>(i);
        const Number angle = Number{2} * static_cast<Number>(kPi) * position
            / static_cast<Number>(kSites);
        state.circulation[i] = forcing + static_cast<Number>(0.01) * std::sin(angle);
        state.humidity[i] = humidityBase + humidityRange
            * (Number{0.5} + Number{0.5} * std::sin(angle));
    }
    return state;
}

template <typename Number>
Rates<Number> evaluateRates(const State<Number>& state) {
    Rates<Number> rates{{
        std::vector<Number>(kSites),
        std::vector<Number>(kSites),
        std::vector<Number>(kSites),
        std::vector<Number>(kSites),
    }};

    const Number forcing = static_cast<Number>(8.0);
    const Number humidityTarget = static_cast<Number>(1.0);
    const Number saturationHumidity = static_cast<Number>(0.88);
    const Number liftThreshold = static_cast<Number>(7.5);
    const Number cloudToRainThreshold = static_cast<Number>(0.001);
    const Number condensationRate = static_cast<Number>(0.08);
    const Number moistureRelaxation = static_cast<Number>(0.12);
    const Number moistureDiffusion = static_cast<Number>(0.035);
    const Number cloudDissipation = static_cast<Number>(0.16);
    const Number rainEfficiency = static_cast<Number>(0.55);

    for (std::size_t i = 0; i < kSites; ++i) {
        const long site = static_cast<long>(i);
        const Number left = state.circulation[wrapIndex(site - 1)];
        const Number twoLeft = state.circulation[wrapIndex(site - 2)];
        const Number right = state.circulation[wrapIndex(site + 1)];

        // Lorenz-96 is a compact model of chaotic atmospheric-scale motion.
        rates.change.circulation[i] = (right - twoLeft) * left
            - state.circulation[i] + forcing;

        const Number humidityLeft = state.humidity[wrapIndex(site - 1)];
        const Number humidityRight = state.humidity[wrapIndex(site + 1)];
        const Number humidityCurvature = humidityLeft - Number{2} * state.humidity[i]
            + humidityRight;
        const Number lift = std::max(Number{0}, state.circulation[i] - liftThreshold);
        const Number condensation = condensationRate
            * std::max(Number{0}, state.humidity[i] - saturationHumidity) * lift;
        const Number rainRate = rainEfficiency
            * std::max(Number{0}, state.cloud[i] - cloudToRainThreshold);

        rates.change.humidity[i] = moistureRelaxation * (humidityTarget - state.humidity[i])
            + moistureDiffusion * humidityCurvature - condensation;
        rates.change.cloud[i] = condensation - cloudDissipation * state.cloud[i] - rainRate;
        rates.change.rain[i] = rainRate;
    }
    return rates;
}

template <typename Number>
State<Number> addScaled(const State<Number>& state, const State<Number>& change, Number scale) {
    State<Number> result = state;
    for (std::size_t i = 0; i < kSites; ++i) {
        result.circulation[i] += scale * change.circulation[i];
        result.humidity[i] += scale * change.humidity[i];
        result.cloud[i] += scale * change.cloud[i];
        result.rain[i] += scale * change.rain[i];
    }
    return result;
}

template <typename Number>
void advanceOneStep(State<Number>& state, Number timeStep) {
    const Rates<Number> first = evaluateRates(state);
    const Rates<Number> second = evaluateRates(addScaled(state, first.change, timeStep / Number{2}));
    const Rates<Number> third = evaluateRates(addScaled(state, second.change, timeStep / Number{2}));
    const Rates<Number> fourth = evaluateRates(addScaled(state, third.change, timeStep));

    for (std::size_t i = 0; i < kSites; ++i) {
        state.circulation[i] += (timeStep / Number{6}) * (first.change.circulation[i]
            + Number{2} * second.change.circulation[i]
            + Number{2} * third.change.circulation[i] + fourth.change.circulation[i]);
        state.humidity[i] += (timeStep / Number{6}) * (first.change.humidity[i]
            + Number{2} * second.change.humidity[i]
            + Number{2} * third.change.humidity[i] + fourth.change.humidity[i]);
        state.cloud[i] += (timeStep / Number{6}) * (first.change.cloud[i]
            + Number{2} * second.change.cloud[i]
            + Number{2} * third.change.cloud[i] + fourth.change.cloud[i]);
        state.rain[i] += (timeStep / Number{6}) * (first.change.rain[i]
            + Number{2} * second.change.rain[i]
            + Number{2} * third.change.rain[i] + fourth.change.rain[i]);

        // Keep small roundoff dips from making these reservoirs negative.
        state.humidity[i] = std::max(Number{0}, state.humidity[i]);
        state.cloud[i] = std::max(Number{0}, state.cloud[i]);
        state.rain[i] = std::max(Number{0}, state.rain[i]);
    }
}

template <typename Number>
RunResult<Number> simulate(State<Number> initial, std::size_t steps, Number timeStep) {
    for (std::size_t step = 0; step < steps; ++step) {
        advanceOneStep(initial, timeStep);
    }

    RunResult<Number> result{initial, initial.rain, std::vector<bool>(kSites, false)};
    const Number stormRainThreshold = static_cast<Number>(0.08);
    for (std::size_t i = 0; i < kSites; ++i) {
        result.stormBySite[i] = result.rainBySite[i] >= stormRainThreshold;
    }
    return result;
}

template <typename First, typename Second>
double finalStateDifference(const State<First>& first, const State<Second>& second) {
    double sum = 0.0;
    for (std::size_t i = 0; i < kSites; ++i) {
        const double circulationDifference = static_cast<double>(first.circulation[i])
            - static_cast<double>(second.circulation[i]);
        const double humidityDifference = static_cast<double>(first.humidity[i])
            - static_cast<double>(second.humidity[i]);
        const double cloudDifference = static_cast<double>(first.cloud[i])
            - static_cast<double>(second.cloud[i]);
        sum += circulationDifference * circulationDifference
            + humidityDifference * humidityDifference
            + cloudDifference * cloudDifference;
    }
    return std::sqrt(sum / (3.0 * static_cast<double>(kSites)));
}

template <typename First, typename Second>
std::pair<std::size_t, double> compareRain(const RunResult<First>& reference,
                                           const RunResult<Second>& comparison) {
    std::size_t changedStormSites = 0;
    double rainfallDifference = 0.0;
    for (std::size_t i = 0; i < kSites; ++i) {
        if (reference.stormBySite[i] != comparison.stormBySite[i]) {
            ++changedStormSites;
        }
        rainfallDifference += std::abs(static_cast<double>(reference.rainBySite[i])
            - static_cast<double>(comparison.rainBySite[i]));
    }
    return {changedStormSites, rainfallDifference};
}

template <typename First, typename Second>
std::string changedStormSiteIds(const RunResult<First>& reference,
                                const RunResult<Second>& comparison) {
    std::string siteIds;
    for (std::size_t i = 0; i < kSites; ++i) {
        if (reference.stormBySite[i] == comparison.stormBySite[i]) {
            continue;
        }
        if (!siteIds.empty()) {
            siteIds += '|';
        }
        siteIds += std::to_string(i);
    }
    return siteIds.empty() ? "none" : siteIds;
}

std::size_t countStorms(const std::vector<bool>& storms) {
    return static_cast<std::size_t>(std::count(storms.begin(), storms.end(), true));
}

template <typename Number>
double totalRain(const RunResult<Number>& result) {
    double total = 0.0;
    for (const Number amount : result.rainBySite) {
        total += static_cast<double>(amount);
    }
    return total;
}

template <typename First, typename Second>
double localStateDifference(const State<First>& first, const State<Second>& second,
                            std::size_t site) {
    const double circulationDifference = static_cast<double>(first.circulation[site])
        - static_cast<double>(second.circulation[site]);
    const double humidityDifference = static_cast<double>(first.humidity[site])
        - static_cast<double>(second.humidity[site]);
    const double cloudDifference = static_cast<double>(first.cloud[site])
        - static_cast<double>(second.cloud[site]);
    const double squaredDifference = circulationDifference * circulationDifference
        + humidityDifference * humidityDifference + cloudDifference * cloudDifference;
    return std::sqrt(squaredDifference / 3.0);
}

template <typename Number>
void writePerturbationRows(std::ofstream& output, const std::string& precision,
                           const RunResult<Number>& baseline, std::size_t steps,
                           Number timeStep, SensitivityAtlas& atlas) {
    for (std::size_t site = 0; site < kSites; ++site) {
        State<Number> perturbedStart = makeInitialState<Number>();
        const Number original = perturbedStart.circulation[site];
        perturbedStart.circulation[site] = std::nextafter(
            original, std::numeric_limits<Number>::infinity());
        const double ulpDelta = static_cast<double>(perturbedStart.circulation[site])
            - static_cast<double>(original);
        const RunResult<Number> perturbed = simulate(std::move(perturbedStart), steps, timeStep);
        const auto [changedStormSites, rainfallDifference] = compareRain(baseline, perturbed);
        for (std::size_t affectedSite = 0; affectedSite < kSites; ++affectedSite) {
            const double difference = localStateDifference(
                baseline.finalState, perturbed.finalState, affectedSite);
            atlas.stateDifference[affectedSite][site] = difference;
            atlas.stormStatusChanged[affectedSite][site] =
                baseline.stormBySite[affectedSite] != perturbed.stormBySite[affectedSite];
            atlas.largestStateDifference = std::max(atlas.largestStateDifference, difference);
        }

        output << "one_ulp_perturbation," << precision << ',' << site << ',' << ulpDelta << ','
               << steps << ','
               << static_cast<double>(timeStep) << ','
               << countStorms(baseline.stormBySite) << ','
               << countStorms(perturbed.stormBySite) << ',' << changedStormSites << ','
               << changedStormSiteIds(baseline, perturbed) << ','
               << totalRain(baseline) << ',' << totalRain(perturbed) << ',' << rainfallDifference << ','
               << finalStateDifference(baseline.finalState, perturbed.finalState) << '\n';
    }
}

std::string cellColor(double value, double maximum) {
    const double fraction = maximum > 0.0 ? std::clamp(value / maximum, 0.0, 1.0) : 0.0;
    const int red = static_cast<int>(std::lround(246.0 - 163.0 * fraction));
    const int green = static_cast<int>(std::lround(246.0 - 185.0 * fraction));
    const int blue = static_cast<int>(std::lround(250.0 - 25.0 * fraction));
    return "rgb(" + std::to_string(red) + "," + std::to_string(green) + ","
        + std::to_string(blue) + ")";
}

void writeAtlasPanel(std::ofstream& output, const SensitivityAtlas& atlas,
                     const std::string& label, int left, int top, double sharedMaximum) {
    constexpr int cellSize = 16;
    constexpr int gridSize = static_cast<int>(kSites) * cellSize;
    output << "<text x=\"" << left << "\" y=\"" << top - 22
           << "\" class=\"panel-title\">" << label << "</text>\n";
    output << "<text x=\"" << left << "\" y=\"" << top + gridSize + 28
           << "\" class=\"axis-title\">nudged site →</text>\n";
    output << "<text x=\"" << left - 35 << "\" y=\"" << top - 7
           << "\" class=\"axis-title\">affected site</text>\n";

    for (std::size_t site = 0; site < kSites; ++site) {
        if (site % 5 == 0 || site == kSites - 1) {
            const int x = left + static_cast<int>(site) * cellSize + cellSize / 2;
            const int y = top + static_cast<int>(site) * cellSize + cellSize / 2 + 4;
            output << "<text x=\"" << x << "\" y=\"" << top + gridSize + 14
                   << "\" class=\"tick\" text-anchor=\"middle\">" << site << "</text>\n";
            output << "<text x=\"" << left - 12 << "\" y=\"" << y
                   << "\" class=\"tick\" text-anchor=\"end\">" << site << "</text>\n";
        }
    }

    for (std::size_t affected = 0; affected < kSites; ++affected) {
        for (std::size_t nudged = 0; nudged < kSites; ++nudged) {
            const int x = left + static_cast<int>(nudged) * cellSize;
            const int y = top + static_cast<int>(affected) * cellSize;
            const bool changedStorm = atlas.stormStatusChanged[affected][nudged];
            output << "<rect x=\"" << x << "\" y=\"" << y << "\" width=\""
                   << cellSize - 1 << "\" height=\"" << cellSize - 1
                   << "\" rx=\"2\" fill=\""
                   << cellColor(atlas.stateDifference[affected][nudged], sharedMaximum)
                   << "\" stroke=\"" << (changedStorm ? "#f05a28" : "#d4d8df")
                   << "\" stroke-width=\"" << (changedStorm ? "2.4" : "0.5") << "\">"
                   << "<title>nudge site " << nudged << "; affected site " << affected
                   << "; local state RMS " << std::setprecision(6)
                   << atlas.stateDifference[affected][nudged]
                   << (changedStorm ? "; storm threshold changed" : "; storm threshold unchanged")
                   << "</title></rect>\n";
        }
    }
    output << "<rect x=\"" << left << "\" y=\"" << top << "\" width=\""
           << gridSize << "\" height=\"" << gridSize
           << "\" fill=\"none\" stroke=\"#6b7280\" stroke-width=\"1\"/>\n";
}

void writeSvg(const std::string& path, const SensitivityAtlas& doubleAtlas,
              const SensitivityAtlas& floatAtlas) {
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("Could not open SVG file: " + path);
    }

    const double sharedMaximum = std::max(doubleAtlas.largestStateDifference,
                                          floatAtlas.largestStateDifference);
    output << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"900\" height=\"540\" "
              "viewBox=\"0 0 900 540\">\n"
           << "<title>Stormline sensitivity atlas</title>\n"
           << "<desc>Each column is a nudged site; each row is an affected site. "
              "Darker cells show larger local final-state differences. Orange borders mark "
              "sites whose storm-threshold status changed.</desc>\n"
           << "<style>text{font-family:Segoe UI,Arial,sans-serif;fill:#202431}"
              ".title{font-size:24px;font-weight:700}.subtitle{font-size:13px;fill:#555d6b}"
              ".panel-title{font-size:17px;font-weight:600}.axis-title{font-size:12px;fill:#404755}"
              ".tick{font-size:10px;fill:#626a78}.legend{font-size:12px;fill:#404755}</style>\n"
           << "<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n"
           << "<text x=\"40\" y=\"38\" class=\"title\">Stormline: numerical influence map</text>\n"
           << "<text x=\"40\" y=\"62\" class=\"subtitle\">Same model, one upward ULP nudge per run. "
              "Color shows final-state change; orange outlines mark storm-threshold flips.</text>\n";

    writeAtlasPanel(output, doubleAtlas, "double precision", 105, 115, sharedMaximum);
    writeAtlasPanel(output, floatAtlas, "float precision", 525, 115, sharedMaximum);
    output << "<text x=\"105\" y=\"487\" class=\"legend\">lighter: smaller local RMS difference</text>\n"
           << "<text x=\"525\" y=\"487\" class=\"legend\">darker: larger local RMS difference "
           << "(shared scale; max " << std::setprecision(6) << sharedMaximum << ")</text>\n"
           << "<text x=\"105\" y=\"510\" class=\"legend\">Orange border: storm-like threshold changed at that site.</text>\n"
           << "</svg>\n";
}

Options parseOptions(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        if (argument == "--steps" && i + 1 < argc) {
            options.steps = std::stoul(argv[++i]);
        } else if (argument == "--dt" && i + 1 < argc) {
            options.timeStep = std::stod(argv[++i]);
        } else if (argument == "--output" && i + 1 < argc) {
            options.outputPath = argv[++i];
        } else if (argument == "--svg" && i + 1 < argc) {
            options.svgPath = argv[++i];
        } else if (argument == "--help") {
            std::cout << "Usage: stormline [--steps N] [--dt time_step] "
                         "[--output csv_path] [--svg svg_path]\n";
            std::exit(EXIT_SUCCESS);
        } else {
            throw std::invalid_argument("Unknown or incomplete argument: " + argument);
        }
    }
    if (options.steps == 0 || !std::isfinite(options.timeStep) || options.timeStep <= 0.0) {
        throw std::invalid_argument("Steps and time step must be positive.");
    }
    return options;
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        const RunResult<double> doubleBaseline = simulate(
            makeInitialState<double>(), options.steps, options.timeStep);
        const RunResult<float> floatBaseline = simulate(
            makeInitialState<float>(), options.steps, static_cast<float>(options.timeStep));

        std::ofstream output(options.outputPath);
        if (!output) {
            throw std::runtime_error("Could not open output file: " + options.outputPath);
        }
        output << "experiment,precision,perturbed_site,ulp_delta,steps,time_step,baseline_storm_count,"
                  "compared_storm_count,changed_storm_site_count,changed_storm_site_ids,baseline_rain_total,"
                  "compared_rain_total,rainfall_absolute_difference,"
                  "final_state_rms_difference\n";
        output << std::setprecision(12);

        const auto [precisionChangedSites, precisionRainDifference] =
            compareRain(doubleBaseline, floatBaseline);
        output << "precision_comparison,double_vs_float,all,n/a," << options.steps << ','
               << options.timeStep << ','
               << countStorms(doubleBaseline.stormBySite) << ','
               << countStorms(floatBaseline.stormBySite) << ',' << precisionChangedSites << ','
               << changedStormSiteIds(doubleBaseline, floatBaseline) << ','
               << totalRain(doubleBaseline) << ',' << totalRain(floatBaseline) << ','
               << precisionRainDifference << ','
               << finalStateDifference(doubleBaseline.finalState, floatBaseline.finalState) << '\n';

        SensitivityAtlas doubleAtlas;
        SensitivityAtlas floatAtlas;
        writePerturbationRows(output, "double", doubleBaseline, options.steps,
                              options.timeStep, doubleAtlas);
        writePerturbationRows(output, "float", floatBaseline, options.steps,
                              static_cast<float>(options.timeStep), floatAtlas);
        output.close();
        writeSvg(options.svgPath, doubleAtlas, floatAtlas);

        std::cout << "Double baseline storm-like sites: "
                  << countStorms(doubleBaseline.stormBySite) << '\n'
                  << "Float baseline storm-like sites: "
                  << countStorms(floatBaseline.stormBySite) << '\n'
                  << "Sensitivity rows written: " << (2 * kSites + 1) << '\n'
                  << "Results: " << options.outputPath << '\n'
                  << "Atlas: " << options.svgPath << '\n'
                  << "Model duration: " << options.steps * options.timeStep
                  << " model-time units\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "stormline: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
