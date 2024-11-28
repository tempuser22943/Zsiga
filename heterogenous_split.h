#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>

struct NodeCapabilities {
    float cpu;      // CPU capability in MHz
    float memory;   // Memory capability in GB
    float bandwidth; // Bandwidth capability in Gbps

    float weighted_sum(float cpu_weight, float memory_weight, float bandwidth_weight) const {
        return cpu * cpu_weight + memory * memory_weight + bandwidth * bandwidth_weight;
    }
};

class HeterogeneousSplitter {
    
public:
    int num_procs;
    int total_nodes;
    std::vector<NodeCapabilities> node_capabilities;
    float total_power;
    std::vector<float> power_ratios;
    std::vector<int> node_distribution;
    std::vector<std::pair<int, int>> node_ranges;
    float cpu_weight;
    float memory_weight;
    float bandwidth_weight;

    HeterogeneousSplitter(const std::vector<NodeCapabilities>& node_capabilities,
                          int total_nodes,
                          float cpu_weight = 1.0, float memory_weight = 0.0, float bandwidth_weight = 0.0)
    : num_procs(node_capabilities.size()), total_nodes(total_nodes), node_capabilities(node_capabilities),
      cpu_weight(cpu_weight), memory_weight(memory_weight), bandwidth_weight(bandwidth_weight) {
        calculateTotalPower();
        calculatePowerRatios();
        calculateNodeDistribution();
        calculateNodeRanges();
    }

    void calculateTotalPower() {
        total_power = 0.0f;
        for (const auto& nc : node_capabilities) {
            total_power += nc.weighted_sum(cpu_weight, memory_weight, bandwidth_weight);
        }
    }

    void calculatePowerRatios() {
        power_ratios.resize(num_procs);
        for (int i = 0; i < num_procs; ++i) {
            power_ratios[i] = node_capabilities[i].weighted_sum(cpu_weight, memory_weight, bandwidth_weight) / total_power;
        }
    }

    void calculateNodeDistribution() {
        node_distribution.resize(num_procs);
        int remaining_nodes = total_nodes;
        int assigned_nodes = 0;

        for (int i = 0; i < num_procs - 1; ++i) {
            node_distribution[i] = std::round(power_ratios[i] * total_nodes);
            assigned_nodes += node_distribution[i];
        }

        // Assign remaining nodes to the last processor
        node_distribution[num_procs - 1] = total_nodes - assigned_nodes;
    }

    void calculateNodeRanges() {
        node_ranges.resize(num_procs);
        int start = 0;
        for (int i = 0; i < num_procs; ++i) {
            int end = start + node_distribution[i] - 1;
            node_ranges[i] = {start, end};
            start = end + 1;
        }
    }

    int get_pid_for_node(int node) const {
        if (node < 0 || node >= total_nodes) {
            throw std::out_of_range("Node index out of range");
        }

        for (int pid = 0; pid < num_procs; ++pid) {
            if (node >= node_ranges[pid].first && node <= node_ranges[pid].second) {
                return pid;
            }
        }

        return num_procs - 1; // Fallback to last processor
    }

    int getNodeCountForProcessor(int processor_id) const {
        if (processor_id < 0 || processor_id >= num_procs) {
            return 0; // Or throw an exception
        }
        return node_distribution[processor_id];
    }

    std::pair<int, int> getNodeRangeForProcessor(int processor_id) const {
        if (processor_id < 0 || processor_id >= num_procs) {
            throw std::out_of_range("Processor ID out of range");
        }
        return node_ranges[processor_id];
    }

    float getPowerRatio(int processor_id) const {
        if (processor_id < 0 || processor_id >= num_procs) {
            throw std::out_of_range("Processor ID out of range");
        }
        return power_ratios[processor_id];
    }
};