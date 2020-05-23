#pragma once

class HyperGraph
{
private:
    /// _numV: number of nodes in the graph.
    uint32_t _numV;
    /// _numE: number of edges in the graph.
    size_t _numE;
    /// _numRRsets: number of RR sets.
    size_t _numRRsets = 0;
    std::vector<bool> _vecVisitBool;
    Nodelist _vecVisitNode;

    size_t _hit = 0;
    double _numSamplesEval = 0;
    double _hyperedgeAvgEval = 0.0;

    bool _isVanilla = false;

    /// Initialization
    void InitHypergraph()
    {
        _numV = (uint32_t)_graph.size();

        for (auto& nbrs : _graph) _numE += nbrs.size();

        _FRsets = FRsets(_numV);
        _vecVisitBool = std::vector<bool>(_numV);
        _vecVisitNode = Nodelist(_numV);
    }

public:
    /// _graph: reverse graph
    const Graph& _graph;
    /// _FRsets: forward cover sets, _FRsets[i] is the node sets that node i can reach
    FRsets _FRsets;
    /// _RRsets: reverse cover sets, _RRsets[i] is the node set that can reach node i
    RRsets _RRsets;

    ProbDist _probDist = WEIGHTS;

    explicit HyperGraph(const Graph& graph) : _graph(graph)
    {
        InitHypergraph();
    }

    /// Set cascade model
    void set_prob_dist(const ProbDist dist)
    {
        _probDist = dist;
    }

    void set_vanilla_sample(const bool isVanilla)
    {
        _isVanilla = isVanilla;
    }

    /// Returns the number of nodes in the graph.
    uint32_t get_nodes() const
    {
        return _numV;
    }

    /// Returns the number of edges in the graph.
    size_t get_edges() const
    {
        return _numE;
    }

    /// Returns the number of RR sets in the graph.
    size_t get_RR_sets_size() const
    {
        return _numRRsets;
    }

    /// Generate a set of n RR sets
    void BuildRRsets(const size_t numSamples)
    {
        void (*func)(const uint32_t uStart, const size_t hyperIdx);

        if (numSamples > SIZE_MAX)
        {
            std::cout << "Error:R too large" << std::endl;
            exit(1);
        }

        const auto prevSize = _numRRsets;
        _numRRsets = _numRRsets > numSamples ? _numRRsets : numSamples;

        if (_isVanilla)
        {
            // std::cout << "Sample RR set by vanilla method" << std::endl;

            for (auto i = prevSize; i < numSamples; i++)
            {
                BuildOneRRset(dsfmt_gv_genrand_uint32_range(_numV), i);
            }

            return ;
        }

        if (_probDist == WC)
        {
            // std::cout << "Sample RR sets in WC model" << std::endl;

            for (auto i = prevSize; i < numSamples; i++)
            {
                BuildOneRRsetWeighted(dsfmt_gv_genrand_uint32_range(_numV), i);
            }

            return ;
        }
        else if (_probDist == UNIFORM)
        {
            // std::cout << "Sample RR sets in uniform model" << std::endl;

            for (auto i = prevSize; i < numSamples; i++)
            {
                BuildOneRRsetConstant(dsfmt_gv_genrand_uint32_range(_numV), i);
            }

            return ;
        }
        else if (_probDist == SKEWED || _probDist == WEIGHTS)
        {
            // std::cout << "Sample RR sets in skewed or weights case" << std::endl;

            for (auto i = prevSize; i < numSamples; i++)
            {
                BuildOneRRsetSkewed(dsfmt_gv_genrand_uint32_range(_numV), i);
            }

            return ;
        }
        else
        {
            for (auto i = prevSize; i < numSamples; i++)
            {
                BuildOneRRset(dsfmt_gv_genrand_uint32_range(_numV), i);;
            }
        }

        return ;
    }


    double EvalHyperedgeAvg()
    {
        return _hyperedgeAvgEval;
    }


    void display_hyperedge_stat()
    {
        size_t total_hyperedges = 0;

        for (size_t i = 0; i < _numRRsets; i++)
        {
            total_hyperedges += _RRsets[i].size();
        }

        double ave_hyperedge_size = 1.0 * total_hyperedges / _numRRsets;
        double var = 0.0;

        for (size_t i = 0; i < _numRRsets; i++)
        {
            double diff = _RRsets[i].size() - ave_hyperedge_size;
            var += (diff * diff);
        }

        var = var / _numRRsets;
        std::cout << "final average RR set size: " << ave_hyperedge_size << ", final variance: " << var << std::endl;
        return ;
    }

    double HyperedgeAvg()
    {
        size_t totalHyperedges = 0;
        double avgSize = 0.0;

        for (size_t i = 0; i < _numRRsets; i++)
        {
            totalHyperedges += _RRsets[i].size();
        }

        avgSize = 1.0 * totalHyperedges / _numRRsets;
        return avgSize;
    }

    double HyperedgeMedian()
    {
        std::vector<int> RRsetSize(_numRRsets);
        for (size_t i = 0; i < _numRRsets; i++)
        {
          RRsetSize[i]  = _RRsets[i].size();
          std::cout << RRsetSize[i] << std::endl;
        }

        std::sort(RRsetSize.begin(), RRsetSize.end());
        int index = _numRRsets/2;
        if (_numRRsets%2 == 1)
        {
            return RRsetSize[index];
        }
        else
        {
            return (RRsetSize[index] + RRsetSize[index-1])/2.0;
        }
    } 

    // Generate one RR set
    void BuildOneRRset(const uint32_t uStart, const size_t hyperIdx)
    {
        size_t numVisitNode = 0, currIdx = 0;
        _FRsets[uStart].push_back(hyperIdx);
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;

        while (currIdx < numVisitNode)
        {
            const auto expand = _vecVisitNode[currIdx++];

            for (auto& nbr : _graph[expand])
            {
                const auto nbrId = nbr.first;

                if (_vecVisitBool[nbrId])
                    continue;

                const auto randDouble = dsfmt_gv_genrand_open_close();

                if (randDouble > nbr.second)
                    continue;

                _vecVisitNode[numVisitNode++] = nbrId;
                _vecVisitBool[nbrId] = true;
                _FRsets[nbrId].push_back(hyperIdx);
            }
        }

        for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;

        _RRsets.push_back(RRset(_vecVisitNode.begin(), _vecVisitNode.begin() + numVisitNode));
    }

    // Independent cascade with weighted probability
    void BuildOneRRsetWeighted(const uint32_t uStart, const size_t hyperIdx)
    {
        size_t numVisitNode = 0, currIdx = 0;
        _FRsets[uStart].push_back(hyperIdx);
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;

        while (currIdx < numVisitNode)
        {
            const auto expand = _vecVisitNode[currIdx++];
            //if (_cascadeModel == IC)
            {
                if (_graph[expand].size() == 0) continue;

                double p =  _graph[expand][0].second;
                double log2Prob = Logarithm(1 - p);

                if (p < 1)
                {
                    double prob = dsfmt_gv_genrand_open_close();
                    int startPos = Logarithm(prob) / log2Prob;
                    int endPos = _graph[expand].size();

                    while (startPos < endPos)
                    {
                        const auto nbrId = _graph[expand][startPos].first;

                        if (_vecVisitBool[nbrId])
                        {
                            int increment = Logarithm(dsfmt_gv_genrand_open_close()) / log2Prob;
                            startPos += (increment + 1);
                            continue;
                        }

                        _vecVisitNode[numVisitNode++] = nbrId;
                        _vecVisitBool[nbrId] = true;
                        _FRsets[nbrId].push_back(hyperIdx);
                        int increment = Logarithm(dsfmt_gv_genrand_open_close()) / log2Prob;
                        startPos += increment + 1;
                    }
                }
                else
                {
                    for (auto& nbr : _graph[expand])
                    {
                        const auto nbrId = nbr.first;

                        if (_vecVisitBool[nbrId])
                            continue;

                        _vecVisitNode[numVisitNode++] = nbrId;
                        _vecVisitBool[nbrId] = true;
                        _FRsets[nbrId].push_back(hyperIdx);
                    }
                }
            }
        }

        for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;

        _RRsets.push_back(RRset(_vecVisitNode.begin(), _vecVisitNode.begin() + numVisitNode));
    }

    /* independent cascade with constant probability */
    void BuildOneRRsetConstant(const uint32_t uStart, const size_t hyperIdx)
    {
        size_t numVisitNode = 0, currIdx = 0;
        _FRsets[uStart].push_back(hyperIdx);
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;
        const double p =  _graph[0][0].second;
        const double const_prob = Logarithm(1 - p);

        if (p == 1)
        {
            while (currIdx < numVisitNode)
            {
                const auto expand = _vecVisitNode[currIdx++];

                //std::cout<<_graph[expand].size()<<std::endl;
                if (_graph[expand].size() == 0) continue;

                for (auto& nbr : _graph[expand])
                {
                    const auto nbrId = nbr.first;

                    //std::cout<<nbr.first<<" "<<nbr.second<<std::endl;
                    if (_vecVisitBool[nbrId])
                        continue;

                    _vecVisitNode[numVisitNode++] = nbrId;
                    _vecVisitBool[nbrId] = true;
                    _FRsets[nbrId].push_back(hyperIdx);
                }
            }
        }
        else
        {
            while (currIdx < numVisitNode)
            {
                const auto expand = _vecVisitNode[currIdx++];

                if (0 == _graph[expand].size())
                {
                    continue;
                }

                int startPos = Logarithm(dsfmt_gv_genrand_open_close()) / const_prob;
                int endPos = _graph[expand].size();

                while (startPos < endPos)
                {
                    //std::cout<<"enter loop"<<std::endl;
                    const auto nbrId = _graph[expand][startPos].first;

                    if (!_vecVisitBool[nbrId])
                    {
                        _vecVisitNode[numVisitNode++] = nbrId;
                        _vecVisitBool[nbrId] = true;
                        _FRsets[nbrId].push_back(hyperIdx);
                    }

                    int increment = Logarithm(dsfmt_gv_genrand_open_close()) / const_prob;
                    startPos += increment + 1;
                }
            }
        }

        for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;

        _RRsets.push_back(RRset(_vecVisitNode.begin(), _vecVisitNode.begin() + numVisitNode));
    }

    // independent cascade with skewed distribution
    void BuildOneRRsetSkewed(const uint32_t uStart, const size_t hyperIdx)
    {
        size_t numVisitNode = 0, currIdx = 0;
        _FRsets[uStart].push_back(hyperIdx);
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;
        double p_threshold = 0.1;

        while (currIdx < numVisitNode)
        {
            const auto expand = _vecVisitNode[currIdx++];
            size_t out_degree = _graph[expand].size();

            if (out_degree > 0)
            {
                size_t startMin = 0;
                size_t endMax = out_degree;

                while (startMin < endMax)
                {
                    const auto &currentedge = _graph[expand][startMin];
                    const auto node_prob = currentedge.second;
                    const auto nbrId = currentedge.first;

                    if (node_prob < p_threshold)
                    {
                        break;
                    }

                    const auto randDouble = dsfmt_gv_genrand_open_close();
                    startMin++;

                    if (randDouble > node_prob) continue;

                    if (_vecVisitBool[nbrId]) continue;

                    _vecVisitNode[numVisitNode++] =  nbrId;
                    _vecVisitBool[nbrId]  = true;
                    _FRsets[nbrId].push_back(hyperIdx);
                }

                while (startMin < endMax)
                {
                    double bucket_probability = _graph[expand][startMin].second;
                    const double log_prob = Logarithm(1 - bucket_probability);
                    double prob = dsfmt_gv_genrand_open_close();
                    startMin += floor(Logarithm(prob) / log_prob);

                    if (startMin >= endMax)
                    {
                        break;
                    }

                    const auto &currentedge = _graph[expand][startMin];
                    const auto nbrId = currentedge.first;
                    const auto accept_probability = currentedge.second;
                    double randDouble = dsfmt_gv_genrand_open_close();
                    startMin++;

                    if (randDouble > accept_probability / bucket_probability || _vecVisitBool[nbrId])
                    {
                        continue;
                    }

                    _vecVisitNode[numVisitNode++] = nbrId;
                    _vecVisitBool[nbrId] = true;
                    _FRsets[nbrId].push_back(hyperIdx);
                }
            }
        }

        for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;

        _RRsets.push_back(RRset(_vecVisitNode.begin(), _vecVisitNode.begin() + numVisitNode));
    }

    // Evaluate the influence spread of a seed set on current generated RR sets
    double CalculateInf(const Nodelist& vecSeed)
    {
        std::vector<bool> vecBoolVst = std::vector<bool>(_numRRsets);
        std::vector<bool> vecBoolSeed(_numV);

        for (auto seed : vecSeed) vecBoolSeed[seed] = true;

        for (auto seed : vecSeed)
        {
            for (auto node : _FRsets[seed])
            {
                vecBoolVst[node] = true;
            }
        }

        size_t count = std::count(vecBoolVst.begin(), vecBoolVst.end(), true);
        return 1.0 * count * _numV / _numRRsets;
    }

    // Efficiently estimate the influence spread with sampling error epsilon within probability 1-delta
    double EfficInfVldtAlg(const Nodelist& vecSeed, const double delta = 1e-3, const double eps = 0.01)
    {
        const double c = 2.0 * (exp(1.0) - 2.0);
        const double LambdaL = 1.0 + 2.0 * c * (1.0 + eps) * log(2.0 / delta) / (eps * eps);
        size_t numHyperEdge = 0;
        size_t numCoverd = 0;
        std::vector<bool> vecBoolSeed(_numV);

        for (auto seed : vecSeed) vecBoolSeed[seed] = true;

        while (numCoverd < LambdaL)
        {
            numHyperEdge++;
            size_t numVisitNode = 0, currIdx = 0;
            const auto uStart = dsfmt_gv_genrand_uint32_range(_numV);

            if (vecBoolSeed[uStart])
            {
                // Stop, this sample is covered
                numCoverd++;
                continue;
            }

            _vecVisitNode[numVisitNode++] = uStart;
            _vecVisitBool[uStart] = true;

            while (currIdx < numVisitNode)
            {
                const auto expand = _vecVisitNode[currIdx++];

                for (auto& nbr : _graph[expand])
                {
                    const auto nbrId = nbr.first;

                    if (_vecVisitBool[nbrId])
                        continue;

                    const auto randDouble = dsfmt_gv_genrand_open_close();

                    if (randDouble > nbr.second)
                        continue;

                    if (vecBoolSeed[nbrId])
                    {
                        // Stop, this sample is covered
                        numCoverd++;
                        goto postProcess;
                    }

                    _vecVisitNode[numVisitNode++] = nbrId;
                    _vecVisitBool[nbrId] = true;
                }
            }

postProcess:

            for (auto i = 0; i < numVisitNode; i++)
                _vecVisitBool[_vecVisitNode[i]] = false;
        }

        return 1.0 * numCoverd * _numV / numHyperEdge;
    }

    // Refresh the hypergraph
    void RefreshHypergraph()
    {
        if (_RRsets.size() != 0)
        {
            for (auto i = _numRRsets; i--;)
            {
                RRset().swap(_RRsets[i]);
            }

            RRsets().swap(_RRsets);

            for (auto i = _numV; i--;)
            {
                FRset().swap(_FRsets[i]);
            }
        }

        _numRRsets = 0;
        _hit = 0;
    }

    // Release memory
    void ReleaseMemory()
    {
        RefreshHypergraph();
        std::vector<bool>().swap(_vecVisitBool);
        Nodelist().swap(_vecVisitNode);
        FRsets().swap(_FRsets);
    }

    void BuildOneRRsetEarlyStopByVanilla(std::unordered_set<uint32_t> &connSet, const uint32_t uStart, const size_t hyperIdx)
    {
        size_t numVisitNode = 0, currIdx = 0;
        _FRsets[uStart].push_back(hyperIdx);
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;

        if (connSet.find(uStart) != connSet.end())
        {
            _hit++;
            goto finished;
        }

        while (currIdx < numVisitNode)
        {
            const auto expand = _vecVisitNode[currIdx++];

            if (0 == _graph[expand].size())
            {
                continue;
            }

            for (auto& nbr : _graph[expand])
            {
                const auto nbrId = nbr.first;

                if (_vecVisitBool[nbrId])
                    continue;

                const auto randDouble = dsfmt_gv_genrand_open_close();

                if (randDouble > nbr.second)
                    continue;

                _vecVisitNode[numVisitNode++] = nbrId;
                _vecVisitBool[nbrId] = true;
                _FRsets[nbrId].push_back(hyperIdx);

                if (connSet.find(nbrId) != connSet.end())
                {
                    _hit++;
                    goto finished;
                }
            }
        }

finished:

        for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;

        _RRsets.push_back(RRset(_vecVisitNode.begin(), _vecVisitNode.begin() + numVisitNode));
    }


    void BuildOneRRsetEarlyStopBySubsim(std::unordered_set<uint32_t> &connSet, const uint32_t uStart, const size_t hyperIdx)
    {
        size_t numVisitNode = 0, currIdx = 0;
        _FRsets[uStart].push_back(hyperIdx);
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;

        if (connSet.find(uStart) != connSet.end())
        {
            _hit++;
            goto finished;
        }

        while (currIdx < numVisitNode)
        {
            const auto expand = _vecVisitNode[currIdx++];
            const double p =  _graph[expand][0].second;

            if (0 == _graph[expand].size())
            {
                continue;
            }

            if (p >= 1.0)
            {
                for (auto& nbr : _graph[expand])
                {
                    const auto nbrId = nbr.first;

                    if (_vecVisitBool[nbrId])
                        continue;

                    _vecVisitNode[numVisitNode++] = nbrId;
                    _vecVisitBool[nbrId] = true;
                    _FRsets[nbrId].push_back(hyperIdx);
                }

                continue;
            }

            const double const_prob = Logarithm(1 - p);
            int startPos = Logarithm(dsfmt_gv_genrand_open_close()) / const_prob;
            int endPos = _graph[expand].size();

            while (startPos < endPos)
            {
                const auto nbrId = _graph[expand][startPos].first;

                if (!_vecVisitBool[nbrId])
                {
                    _vecVisitNode[numVisitNode++] = nbrId;
                    _vecVisitBool[nbrId] = true;
                    _FRsets[nbrId].push_back(hyperIdx);
                }

                if (connSet.find(nbrId) != connSet.end())
                {
                    _hit++;
                    goto finished;
                }

                int increment = Logarithm(dsfmt_gv_genrand_open_close()) / const_prob;
                startPos += increment + 1;
            }
        }

finished:

        for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;

        _RRsets.push_back(RRset(_vecVisitNode.begin(), _vecVisitNode.begin() + numVisitNode));
    }

    void BuildRRsetsEarlyStop(std::unordered_set<uint32_t> &connSet, const int numSamples)
    {
        const auto prevSize = _numRRsets;
        _numRRsets = _numRRsets > numSamples ? _numRRsets : numSamples;

        if (_isVanilla)
        {
            // std::cout << "Sample RR sets with early stop by Vanilla method" << std::endl;
            for (auto i = prevSize; i < numSamples; i++)
            {
                BuildOneRRsetEarlyStopByVanilla(connSet, dsfmt_gv_genrand_uint32_range(_numV), i);
            }
        }
        else
        {
            // std::cout << "Sample RR sets with early stop By SUBSIM" << std::endl;
            for (auto i = prevSize; i < numSamples; i++)
            {
                BuildOneRRsetEarlyStopBySubsim(connSet, dsfmt_gv_genrand_uint32_range(_numV), i);
            }
        }
    }

    double EvalSeedSetInfByVanilla(std::unordered_set<uint32_t> &connSet, const int numSamples)
    {
        uint32_t numCovered = 0;
        int64_t totalHyperedgeSize = 0;
        _numSamplesEval = numSamples;

        for (int i = 1; i < numSamples; i++)
        {
            uint32_t uStart = dsfmt_gv_genrand_uint32_range(_numV);
            size_t numVisitNode = 0, currIdx = 0;
            _vecVisitNode[numVisitNode++] = uStart;
            _vecVisitBool[uStart] = true;

            if (connSet.find(uStart) != connSet.end())
            {
                numCovered++;
                goto finished;
            }

            while (currIdx < numVisitNode)
            {
                const auto expand = _vecVisitNode[currIdx++];

                if (0 == _graph[expand].size())
                {
                    continue;
                }

                const double p =  _graph[expand][0].second;

                for (auto& nbr : _graph[expand])
                {
                    const auto nbrId = nbr.first;

                    if (_vecVisitBool[nbrId])
                        continue;

                    const auto randDouble = dsfmt_gv_genrand_open_close();

                    if (randDouble > nbr.second)
                        continue;

                    _vecVisitNode[numVisitNode++] = nbrId;
                    _vecVisitBool[nbrId] = true;

                    if (connSet.find(nbrId) != connSet.end())
                    {
                        numCovered++;
                        goto finished;
                    }
                }
            }

finished:
            totalHyperedgeSize += numVisitNode;

            for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;
        }

        _hyperedgeAvgEval = 1.0 * totalHyperedgeSize / numSamples;
        return 1.0 * numCovered * _numV / numSamples;
    }

    double EvalSeedSetInfBySubsim(std::unordered_set<uint32_t> &connSet, const int numSamples)
    {
        uint32_t numCovered = 0;
        int64_t totalHyperedgeSize = 0;
        _numSamplesEval = numSamples;

        for (int i = 1; i < numSamples; i++)
        {
            uint32_t uStart = dsfmt_gv_genrand_uint32_range(_numV);
            size_t numVisitNode = 0, currIdx = 0;
            _vecVisitNode[numVisitNode++] = uStart;
            _vecVisitBool[uStart] = true;

            if (connSet.find(uStart) != connSet.end())
            {
                numCovered++;
                goto finished;
            }

            while (currIdx < numVisitNode)
            {
                const auto expand = _vecVisitNode[currIdx++];

                if (0 == _graph[expand].size())
                {
                    continue;
                }

                const double p =  _graph[expand][0].second;

                if (p >= 1.0)
                {
                    for (auto& nbr : _graph[expand])
                    {
                        const auto nbrId = nbr.first;

                        if (_vecVisitBool[nbrId])
                            continue;

                        _vecVisitNode[numVisitNode++] = nbrId;
                        _vecVisitBool[nbrId] = true;

                        if (connSet.find(nbrId) != connSet.end())
                        {
                            numCovered++;
                            goto finished;
                        }
                    }

                    continue;
                }

                const double const_prob = Logarithm(1 - p);
                int startPos = Logarithm(dsfmt_gv_genrand_open_close()) / const_prob;
                int endPos = _graph[expand].size();

                while (startPos < endPos)
                {
                    const auto nbrId = _graph[expand][startPos].first;

                    if (!_vecVisitBool[nbrId])
                    {
                        _vecVisitNode[numVisitNode++] = nbrId;
                        _vecVisitBool[nbrId] = true;
                    }

                    if (connSet.find(nbrId) != connSet.end())
                    {
                        numCovered++;
                        goto finished;
                    }

                    int increment = Logarithm(dsfmt_gv_genrand_open_close()) / const_prob;
                    startPos += increment + 1;
                }
            }

finished:
            totalHyperedgeSize += numVisitNode;

            for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;
        }

        _hyperedgeAvgEval = 1.0 * totalHyperedgeSize / numSamples;
        return 1.0 * numCovered * _numV / numSamples;
    }

    double EvalSeedSetInf(std::unordered_set<uint32_t> &connSet, const int numSamples)
    {
        if (_isVanilla)
        {
            return EvalSeedSetInfByVanilla(connSet, numSamples);
        }
        else
        {
            return EvalSeedSetInfBySubsim(connSet, numSamples);
        }
    }

    double CalculateInfEarlyStop()
    {
        return 1.0 * _hit * _numV / _numRRsets;
    }
};

using THyperGraph = HyperGraph;
using PHyperGraph = std::shared_ptr<THyperGraph>;
