#include "heterogenous_hybridmap.h"
#include "MessageManager.h"
#include <unordered_map>
#include <vector>
#include <iostream>
#include <algorithm>
#include "ska/flat_hash_map.hpp"

#define pure(x) ((x >> 31) ^ x)

using std::swap;

namespace zsiga {

void _part(Edge a[], int l, int r, int &i, int &j) {
    if (r - l <= 1) {
        if (a[r].v < a[l].v)
            swap(a[r], a[l]);
        i = l;
        j = r;
        return;
    }
    
    int m = l;
    int pivot = a[r].v;
    while (m <= r) {
        if (a[m].v < pivot)
            swap(a[l++], a[m++]);
        else if (a[m].v > pivot)
            swap(a[m], a[r--]);
        else
            m++;
    }
    i = l-1;
    j = m;
}

void _sort(Edge a[], int l, int r) {
    if (l >= r) return;
    
    int i, j;
    _part(a, l, r, i, j);
    _sort(a, l, i);
    _sort(a, j, r);
}

class RemInit {
public:
    int num_nodes, num_procs;
    hash_t p;
    HeterogeneousSplitter* splitter;

    RemInit(int num_nodes, int num_procs, HeterogeneousSplitter* splitter)
        : num_nodes(num_nodes), num_procs(num_procs), splitter(splitter) {}

    ~RemInit() {
        p.clear();
    }

    int pp(int u) {
        auto x = p.find(u);
        return x == p.end() ? u : x->second;
    }

    void merge(int u, int v) {
        int up, vp;
        while ((up = pp(u)) != (vp = pp(v))) {
            if (up < vp) {
                p[v] = up;
                v = vp;
            } else {
                p[u] = vp;
                v = up;
            }
        }
    }

    int find(int u) {
        int up = pp(u);
        int upp;
        while (up != (upp = pp(up))) {
            p[u] = upp;
            u = up;
            up = upp;
        }
        return up;
    }

    int partition(Edge edges_tmp[]) {
        int size = 0;
        for (auto pair : p) {
            edges_tmp[size++] = {pair.first, find(pair.second)};
        }

        _sort(edges_tmp, 0, size-1);
        
        std::vector<int> heads(num_procs, -1);
        
        for (int i = 0, j = 0; i < size;) {
            std::fill(heads.begin(), heads.end(), -1);
            int cc = edges_tmp[i].v;
            int cc_pid = splitter->get_pid_for_node(cc);
            heads[cc_pid] = cc;
            for (; j < size && edges_tmp[j].v == cc; j++) {
                int jp = splitter->get_pid_for_node(edges_tmp[j].u);
                if (heads[jp] > edges_tmp[j].u || heads[jp] == -1) heads[jp] = edges_tmp[j].u;
            }

            for (; i < j; i++) {
                int ip = splitter->get_pid_for_node(edges_tmp[i].u);
                if (heads[ip] != edges_tmp[i].u)
                    edges_tmp[i].v = heads[ip];
            }
        }

        return size; 
    }

    void clear() {
        p.clear();
    }
};

class RemZsiga {
public:
    HybridMap p;
    int num_procs;
    HeterogeneousSplitter* splitter;

    RemZsiga(int num_nodes, int pid, HeterogeneousSplitter* splitter)
        : p(num_nodes, pid, splitter), num_procs(splitter->num_procs), splitter(splitter) {}

    void merge(int u, int v) {
        int up = pure(p[u]);
        int vp = pure(p[v]);

        if (up == u && vp == v) {
            p[u] = vp;
            return;
        }

        while (up != vp) {
            if (up < vp) {
                p[v] = ~up;
                v = vp;
            } else {
                p[u] = ~vp;
                u = up;
            }
            up = pure(p[u]);
            vp = pure(p[v]);
        }
    }

    int find(int u) {
        int up = pure(p[u]);
        int upp;
        while (up != (upp = pure(p[up]))) {
            p[u] = ~upp;
            u = up;
            up = upp;
        }
        return up;
    }

    int getEdgesToEmit(Edge toemit[], size_t& num_changes) {
        int size = 0;

        for (const auto &x : p.p_out) {
            int pr = find(x.first);
            if (pr != x.first && x.second < 0) {
                toemit[size++] = Edge{x.first, pr};
                num_changes++;
            }
        }

        _sort(toemit, 0, size-1);
        
        std::vector<int> heads(num_procs, -1);
        
        for (int i = 0, j = 0; i < size;) {
            std::fill(heads.begin(), heads.end(), -1);
            int cc = toemit[i].v;
            int cc_pid = splitter->get_pid_for_node(cc);
            heads[cc_pid] = cc;
            for (; j < size && toemit[j].v == cc; ++j) {
                int jp = splitter->get_pid_for_node(toemit[j].u);
                if (heads[jp] > toemit[j].u || heads[jp] == -1) heads[jp] = toemit[j].u;
            }

            for (; i < j; ++i) {
                int ip = splitter->get_pid_for_node(toemit[i].u);
                if (heads[ip] != toemit[i].u)
                    toemit[i].v = heads[ip];
            }
        }
        

        ska::flat_hash_map<int, int> ucc;

        std::pair<int, int> node_range = splitter->getNodeRangeForProcessor(p.pid);

        for (int node=node_range.first; node<=node_range.second; node++) {
            // int u = p.p_in[i]; //parent of i th node
            int up = find(node); //global root
            if (splitter->get_pid_for_node(up) == p.pid) continue; //global root is the local root
            auto it = ucc.find(up);
            if (it == ucc.end()) {
                ucc[up] = node;
                toemit[size++] = Edge{node, up};
                if (p[node] < 0) {
                    p[node] = up;
                    num_changes++;
                }
            } else {
                p[node] = it->second;
            }
        }

        p.p_out.clear();
        
        return size;
    }

    void _sort(Edge a[], int l, int r){
        if (l >= r) return;
        
            int i, j;
            _part(a, l, r, i, j);
            _sort(a, l, i);
            _sort(a, j, r);
        }

        void _part(Edge a[], int l, int r, int &i, int &j){
            if (r - l <= 1) {
                if (a[r].v < a[l].v)
                    swap(a[r], a[l]);
                i = l;
                j = r;
                return;
            }
            
            int m = l;
            int pivot = a[r].v;
            while (m <= r) {
                if (a[m].v < pivot)
                    swap(a[l++], a[m++]);
                else if (a[m].v > pivot)
                    swap(a[m], a[r--]);
                else
                    m++;
            }
            i = l-1;
            j = m;
        }

    };
}