#pragma once


class Argument
{
public:
    // Function parameter.
    // format: format graph
    // im: influence maximization
    std::string _funcStr = "im";
    FuncType _func = IM;

    // The number of nodes to be selected. Default is 50.
    int _seedsize = 50;

    // For the uniform setting, every edge has the same diffusion probability.
    float _probEdge = float(0.1);

    // Error threshold 1-1/e-epsilon.
    double _eps = 0.1;

    // Failure probability delta. Default is 1/#nodes.
    double _delta = -1.0;

    // Graph name. Default is "facebook".
    std::string _graphname = "facebook";

    // Probability distribution
    // weights: graph data with weights
    // wc: wc setting
    // uniform: uniform setting
    // skewed: skewed distribution
    std::string _probDistStr = "wc";
    ProbDist _probDist = WC;
    std::string _skewType = "exp";

    // Directory
    std::string _dir = "graphInfo";

    // Result folder
    std::string _resultFolder = "result";

    // File name of the result
    std::string _outFileName;

    // wc variant
    double _wcVar = 1.0;

    // sample RR set with the vanilla method
    bool _vanilla = false;

    // use hist algorithm
    bool _hist = false;

    Argument(int argc, char* argv[])
    {
        std::string param, value;

        for (int ind = 1; ind < argc; ind++)
        {
            if (argv[ind][0] != '-') break;

            std::stringstream sstr(argv[ind]);
            getline(sstr, param, '=');
            getline(sstr, value, '=');

            if (!param.compare("-func")) _funcStr = value;
            else if (!param.compare("-seedsize")) _seedsize = stoi(value);
            else if (!param.compare("-eps")) _eps = stod(value);
            else if (!param.compare("-delta")) _delta = stod(value);
            else if (!param.compare("-gname")) _graphname = value;
            else if (!param.compare("-dir")) _dir = value;
            else if (!param.compare("-outpath")) _resultFolder = value;
            else if (!param.compare("-pdist")) _probDistStr = value;
            else if (!param.compare("-pedge")) _probEdge = stof(value);
            else if (!param.compare("-wcvariant")) _wcVar = stod(value);
            else if (!param.compare("-skew")) _skewType = value;
            else if (!param.compare("-vanilla")) _vanilla = (value == "1");
            else if (!param.compare("-hist")) _hist = (value == "1");
        }

        if (_wcVar <= 0)
        {
            //wrong input
            _wcVar = 1.0;
        }

        decode_func_type();
        decode_prob_dist();
    }



    void build_outfilename(int seedSize, ProbDist dist, Graph& graph)
    {
        std::string distStr; 

        if (dist == WEIGHTS)
        {
            _probDistStr = "weights";
        }
        else if (dist == WC)
        {
            _probDistStr = "wc";
        }
        else if (dist == UNIFORM)
        {
            _probDistStr = "uniform";

            for (int i = 0; i < graph.size(); i++)
            {
                if (graph[i].size() > 0 )
                {
                    _probEdge = graph[i][0].second;
                    break;
                }
            }
        }
        else
        {
            _probDistStr = "skewed";
        }

        _outFileName = TIO::BuildOutFileName(_graphname, "subsim", seedSize, _probDistStr, _probEdge);

        return ;
    }

    void decode_prob_dist()
    {
        if (_probDistStr == "wc")
        {
            _probDist = WC;
        }
        else if (_probDistStr == "uniform")
        {
            _probDist = UNIFORM;
        }
        else if (_probDistStr == "skewed")
        {
            _probDist = SKEWED;
        }
        else if (_probDistStr == "weights")
        {
            _probDist = WEIGHTS;
        }
        else 
        {
            _probDist = PROB_DIST_ERROR;
        }
    }

    void decode_func_type()
    {
        if (_funcStr == "format")
        {
            _func = FORMAT;
        }
        else if (_funcStr == "im")
        {
            _func = IM;
        }
        else
        {
            _func = FUNC_ERROR;
        }
    }
};

using TArgument = Argument;
using PArgument = std::shared_ptr<TArgument>;
