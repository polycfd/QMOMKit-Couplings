#include "settings.hpp"

#include <stdexcept>

namespace hyqmom_transport
{

Settings read_settings(const int argc, char* argv[])
{
    Settings settings;
    for (int i = 1; i < argc; i += 2)
    {
        const std::string option = argv[i];
        if (i + 1 >= argc)
        {
            throw std::invalid_argument("The option " + option + " needs a value.");
        }
        const std::string value = argv[i + 1];
        if (option == "--case")
        {
            settings.case_name = value;
        }
        else if (option == "--nodes")
        {
            settings.node_count = std::stoul(value);
        }
        else if (option == "--flux")
        {
            if (value == "nodes")
            {
                settings.flux = FluxType::nodes;
            }
            else if (value == "speeds")
            {
                settings.flux = FluxType::speeds;
            }
            else if (value == "hll")
            {
                settings.flux = FluxType::hll;
            }
            else
            {
                throw std::invalid_argument("The flux must be nodes, speeds or hll.");
            }
        }
        else if (option == "--hll-speeds")
        {
            if (value != "local" && value != "global")
            {
                throw std::invalid_argument("The HLL speeds must be local or global.");
            }
            settings.global_hll_speeds = value == "global";
        }
        else if (option == "--cells")
        {
            settings.cells = std::stoul(value);
        }
        else if (option == "--courant")
        {
            settings.courant = std::stod(value);
        }
        else if (option == "--velocity-width")
        {
            settings.velocity_width = std::stod(value);
        }
        else if (option == "--output")
        {
            settings.output = value;
        }
        else
        {
            throw std::invalid_argument("Unknown option " + option + ".");
        }
    }
    if (settings.case_name.empty() || settings.node_count == 0)
    {
        throw std::invalid_argument("Usage: hyqmomTransport1D --case <case> --nodes <N> --flux nodes|speeds|hll "
                                    "[--hll-speeds local|global] [--cells <n>] [--courant <c>] "
                                    "[--velocity-width <s>] [--output <dir>]");
    }
    return settings;
}

std::string flux_label(const Settings& settings)
{
    switch (settings.flux)
    {
    case FluxType::nodes:
        return "nodes";
    case FluxType::speeds:
        return "speeds";
    case FluxType::hll:
        return settings.global_hll_speeds ? "hll-global" : "hll-local";
    }
    return "";
}

} // namespace hyqmom_transport
