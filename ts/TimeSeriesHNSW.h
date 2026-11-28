#pragma once
#include "../hnswlib/hnswlib.h"
#include "TimeSeriesSpace.h"
#include <iostream>
#include <stdexcept>
#include <stdlib.h>
#include <memory>

/*
TimeSeriesHNSW sets up the memory layout for the index by allocating the contiguous Level 0 buffer 
(where all nodes and their 8-byte window indices live), computing the SHG compression level sizes, 
and initializing the pointer array for upper-level connections. It also declares searchKnn, saveIndex, 
and addDataPoint method as the foundation for searching, saving, and inserting nodes.
*/

template<typename dist_t>
class TimeSeriesHNSW : public hnswlib::AlgorithmInterface<dist_t>{
public:
    char *level0_memory_{nullptr};
    size_t size_data_per_element_{0};
    size_t data_rep_size_{0};
    size_t offsetData_;
    unsigned int maxFixLevel_{0}; // base layer in HNSW is layer0
    unsigned int k_{4}; // for determining amount of layers
    TimeSeriesSpace *ts_space_;
    size_t window_size_;
    char *data_rep_memory_{nullptr};
    size_t size_links_per_element_{0};
    int maxlevel_{-1};
    int enterpoint_node_{-1};
    size_t max_elements_{0};
    using tableint = unsigned int;

    char **linkLists_{nullptr};

    std::default_random_engine level_generator_;
    double mult_{0.0};

    TimeSeriesHNSW(
        TimeSeriesSpace *space, 
        size_t window_size, 
        size_t num_windows, 
        size_t M = 16,
        size_t random_seed = 100) {
        
        ts_space_ = space;

        level_generator_.seed(random_seed);
        mult_ = 1.0 / log(1.0 * M);

        size_t M0 = M * 2;          // Anzahl der Kanten ist für das Basislevel (Level 0) am höchsten

        data_rep_size_ = 0;
        maxFixLevel_ = 0;
        int dim = window_size;
        while(dim/k_>=k_){
            maxFixLevel_++;
            dim = (int) ceil(float(dim)/float(k_));
            data_rep_size_+=dim;
        }
        
        
        size_t link_list_size = sizeof(unsigned int) + M0 * sizeof(int);
        offsetData_ = link_list_size;
        size_data_per_element_ = sizeof(size_t) + link_list_size;

        level0_memory_ = (char *) malloc(num_windows * size_data_per_element_);

        if (level0_memory_ == nullptr) {
                throw std::runtime_error("Not enough memory for Level 0");
        }
        
        linkLists_ = (char **) malloc(num_windows * sizeof(void *));
        if (linkLists_ == nullptr){
            free(level0_memory_);
            throw std::runtime_error("Not enough memory for link lists");
        }

        for (size_t i = 0; i < num_windows; i++){
            linkLists_[i] = nullptr;
        }

        window_size_ = window_size;

        max_elements_ = num_windows;

        size_links_per_element_ = sizeof(unsigned int) + M * sizeof(tableint);

        if (data_rep_size_ > 0) {
            data_rep_memory_ = (char *) malloc(num_windows * data_rep_size_ * sizeof(float));
            if (data_rep_memory_ == nullptr) {
                free(level0_memory_);
                free(linkLists_);
                throw std::runtime_error("Not enough memory for SHG representations");
            }
        }
        
        std::cout << "Level 0 memory allocated at: " << (void*) level0_memory_ << " | Node size: " << size_data_per_element_ << " bytes." << std::endl;
    }

    ~TimeSeriesHNSW(){
        if (level0_memory_ != nullptr) free(level0_memory_);
        if (linkLists_ != nullptr) {
            for (size_t i = 0; i < max_elements_; i++) {
                if (linkLists_[i] != nullptr) free(linkLists_[i]);
            }
            free(linkLists_);
        }
        if (data_rep_memory_ != nullptr) free(data_rep_memory_);
    }

    inline char* getDataByInternalId(size_t internal_id) const {
        return level0_memory_ + internal_id * size_data_per_element_ + offsetData_;
    }

    size_t getDatasizeLevel(int level) const{
        if(level == 0) return window_size_;
        else{
            int size_level = window_size_;
            for(int i =1; i<=level; i++){
                size_level = (int) ceil(float(size_level)/float(k_));
            }
            return size_level;
        }
    }

    std::priority_queue<std::pair<dist_t, hnswlib::labeltype>> searchKnn(
        const void* query, 
        unsigned int k,   // FIX: Hier stand vorher size_t k
        size_t ef, 
        hnswlib::BaseFilterFunctor* isIdAllowed = nullptr
    ) const override {
        // Vorerst leerer Dummy-Return, damit es kompiliert
        throw std::runtime_error("searchKnn ist noch nicht implementiert!");
    }

    // 2. Die Speicherfunktion
    void saveIndex(const std::string &location) override {
        throw std::runtime_error("saveIndex ist noch nicht implementiert!");
    }

    int getRandomLevel(double reverse_size) {
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        double r = -log(distribution(level_generator_)) * reverse_size;
        return (int) r;
    }

    unsigned short int getListCount(linklistsizeint * ptr) const {
        return *((unsigned short int *)ptr);
    }

    linklistsizeint *get_linklist(tableint internal_id, int level) const {
        return (linklistsizeint *) (linkLists_[internal_id] + (level - 1) * size_links_per_element_);
    }

    char *getDataRepreByLevel(tableint internalId, int level) const{
        if(level == 0) return getDataByInternalId(internalId);
        size_t size_sum_previous = 0;
        
        for(int i=1; i<level; i++) size_sum_previous+= getDatasizeLevel(i);
        return (data_rep_memory_ + internalId * data_rep_size_*sizeof(float) + size_sum_previous*sizeof(float));
    }

    dist_t getDisByLevel(tableint data1, tableint data2, int level) const {
        
        if(level >= maxFixLevel_){
            float *pVect1 = (float *) getDataRepreByLevel(data1, maxFixLevel_);
            float *pVect2 = (float *) getDataRepreByLevel(data2, maxFixLevel_);
            float t = *pVect1 - *pVect2;
            return t*t;
        }

        float *pVect1 = (float *) getDataRepreByLevel(data1, level);
        float *pVect2 = (float *) getDataRepreByLevel(data2, level);

        int rep_size = getDatasizeLevel(level);
        float res = 0;
        for (int i = 0; i < rep_size; i++) {
            float t = *pVect1 - *pVect2;
            pVect1++;
            pVect2++;
            res += t * t;
        }
        return res;
    }

    tableint mutuallyConnectNewElement(
        const void *data_point,
        tableint cur_c,
        std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> &top_candidates,
        int level,
        bool isUpdate) {
        size_t Mcurmax = level ? maxM_ : maxM0_;

        getNeighborsByHeuristic2(top_candidates, M_, level);

        if (top_candidates.size() > M_)
            throw std::runtime_error("Should be not be more than M_ candidates returned by the heuristic");

        std::vector<tableint> selectedNeighbors;
        selectedNeighbors.reserve(M_);
        while (top_candidates.size() > 0) {
            selectedNeighbors.push_back(top_candidates.top().second);
            top_candidates.pop();
        }

        tableint next_closest_entry_point = selectedNeighbors.back();

        {
            // lock only during the update
            linklistsizeint *ll_cur;
            if (level == 0)
                ll_cur = get_linklist0(cur_c);
            else
                ll_cur = get_linklist(cur_c, level);

            if (*ll_cur && !isUpdate) {
                throw std::runtime_error("The newly inserted element should have blank link list");
            }
            setListCount(ll_cur, selectedNeighbors.size());
            tableint *data = (tableint *) (ll_cur + 1);
            for (size_t idx = 0; idx < selectedNeighbors.size(); idx++) {
                if (data[idx] && !isUpdate)
                    throw std::runtime_error("Possible memory corruption");
                if (level > element_levels_[selectedNeighbors[idx]])
                    throw std::runtime_error("Trying to make a link on a non-existent level");

                data[idx] = selectedNeighbors[idx];
            }
        }

        for (size_t idx = 0; idx < selectedNeighbors.size(); idx++) {

            linklistsizeint *ll_other;
            if (level == 0)
                ll_other = get_linklist0(selectedNeighbors[idx]);
            else
                ll_other = get_linklist(selectedNeighbors[idx], level);

            size_t sz_link_list_other = getListCount(ll_other);

            if (sz_link_list_other > Mcurmax)
                throw std::runtime_error("Bad value of sz_link_list_other");
            if (selectedNeighbors[idx] == cur_c)
                throw std::runtime_error("Trying to connect an element to itself");
            if (level > element_levels_[selectedNeighbors[idx]])
                throw std::runtime_error("Trying to make a link on a non-existent level");

            tableint *data = (tableint *) (ll_other + 1);

            bool is_cur_c_present = false;
            if (isUpdate) {
                for (size_t j = 0; j < sz_link_list_other; j++) {
                    if (data[j] == cur_c) {
                        is_cur_c_present = true;
                        break;
                    }
                }
            }

            // If cur_c is already present in the neighboring connections of `selectedNeighbors[idx]` then no need to modify any connections or run the heuristics.
            if (!is_cur_c_present) {
                if (sz_link_list_other < Mcurmax) {
                    data[sz_link_list_other] = cur_c;
                    setListCount(ll_other, sz_link_list_other + 1);
                } else {
                    // finding the "weakest" element to replace it with the new one

                    //dist_t d_max = getDisByLevel(getDataByInternalId(cur_c), getDataByInternalId(selectedNeighbors[idx]), level);
                    dist_t d_max = getDisByLevel(cur_c, selectedNeighbors[idx], 0);

                    // Heuristic:
                    std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> candidates;
                    candidates.emplace(d_max, cur_c);

                    for (size_t j = 0; j < sz_link_list_other; j++) {
                        
                        // dist_t d_ = getDisByLevel(getDataByInternalId(data[j]), getDataByInternalId(selectedNeighbors[idx])
                        //                         , level);
                        dist_t d_ = getDisByLevel(data[j], selectedNeighbors[idx]
                                                , 0);

                        candidates.emplace(d_, data[j]);
                        // candidates.emplace(
                        //         fstdistfunc_(getDataByInternalId(data[j]), getDataByInternalId(selectedNeighbors[idx]),
                        //                         dist_func_param_), data[j]);
                    }

                    getNeighborsByHeuristic2(candidates, Mcurmax,level);

                    int indx = 0;
                    while (candidates.size() > 0) {
                        data[indx] = candidates.top().second;
                        candidates.pop();
                        indx++;
                    }

                    setListCount(ll_other, indx);
                    // Nearest K:
                    /*int indx = -1;
                    for (int j = 0; j < sz_link_list_other; j++) {
                        dist_t d = fstdistfunc_(getDataByInternalId(data[j]), getDataByInternalId(rez[idx]), dist_func_param_);
                        if (d > d_max) {
                            indx = j;
                            d_max = d;
                        }
                    }
                    if (indx >= 0) {
                        data[indx] = cur_c;
                    } */
                }
            }
        }

        return next_closest_entry_point;
    }
    
    void addDataPoint(size_t window_idx) {

        int curlevel = getRandomLevel(mult_);

        memset(level0_memory_ + window_idx * size_data_per_element_, 0, size_data_per_element_);
        memcpy(getDataByInternalId(window_idx), &window_idx, sizeof(size_t));

        std::vector<float> rep;
        rep.reserve(data_rep_size_ + window_size_);

        const float* raw_window = ts_space_->getWindowPtr(window_idx);
        float mu = ts_space_->getMean(window_idx);
        float sigma = ts_space_->getStd(window_idx);

        float inv_sigma = (sigma > 1e-6f) ? (1.0f / sigma) : 0.0f;
        for (size_t i = 0; i < window_size_; ++i) {
            rep.push_back((raw_window[i] - mu) * inv_sigma);
        }

        int previous_level_pos = 0;

        //TODO: understand wtf is going on here
        for(int cur_lev = 1; cur_lev <= maxFixLevel_; cur_lev++){
            int previous_level_size = getDatasizeLevel(cur_lev-1);

            for(int i = 0; i< previous_level_size; i+=k_){
                float sum = 0;

                if(i+k_ > previous_level_size){
                    //std::cout  << "i+k_ > previous_level_size " << std::endl;
                    for(int j = i; j< previous_level_size; j++) sum += rep[previous_level_pos + j];
                    float var = float(sum) /float(previous_level_size-i);
                    rep.push_back(var);
                }else{
                    //std::cout  << "else i+k_ > previous_level_size " << std::endl;
                    for(int j = i; j<i+k_;j++){
                        sum += rep[previous_level_pos + j];
                    }                    
                    
                    float var = float(sum) /float(k_);
                    rep.push_back(var);
                }

            }

            previous_level_pos+= previous_level_size;
        }

        rep.erase(rep.begin(), rep.begin() + getDatasizeLevel(0));
        
        if (data_rep_size_ > 0) {
            memcpy(data_rep_memory_ + window_idx * data_rep_size_ * sizeof(float), rep.data(), rep.size() * sizeof(float));
        }
        
        if (curlevel) {
            linkLists_[window_idx] = (char *) malloc(size_links_per_element_ * curlevel + 1);
            if (linkLists_[window_idx] == nullptr)
                throw std::runtime_error("Not enough memory: addPoint failed to allocate linklist");
            memset(linkLists_[window_idx], 0, size_links_per_element_ * curlevel + 1);
        }

        int maxlevelcopy = maxlevel_;
        int currObj = enterpoint_node_;

        if ((signed)currObj != -1) {
            if (curlevel < maxlevelcopy) {

                dist_t curdist = getDisByLevel(window_idx, currObj, 0);

                for (int level = maxlevelcopy; level > curlevel; level--) {
                    
                    //curdist = getDisByLevel(getDataRepreByLevel(cur_c, level),getDataRepreByLevel(currObj, level), level);
                    //curdist = getDisByLevel(cur_c, currObj, level);
                    
                    bool changed = true;
                    while (changed) {
                        changed = false;
                        unsigned int *data;
                        //std::unique_lock <std::mutex> lock(link_list_locks_[currObj]);
                        //std::cout << "get_linklist: " << currObj << " " << level << " " << get_linklist(currObj, level) << " " << *get_linklist(currObj, level) <<  " " << getListCount(get_linklist(currObj, level)) <<   std::endl;
                        data = get_linklist(currObj, level);
                        int size = getListCount(data);

                        tableint *datal = (tableint *) (data + 1);
                        //std::cout << "datal: " << *datal << std::endl;
                        for (int i = 0; i < size; i++) {
                            tableint cand = datal[i];
                            if (cand < 0 || cand > max_elements_)
                                std::cout << cand << " cand error";
                                //throw std::runtime_error("cand error");
                            dist_t d = getDisByLevel(window_idx, cand, 0);

                            if (d < curdist) {
                                curdist = d;
                                currObj = cand;
                                changed = true;
                            }
                        }
                    }
                }
            }

            for (int level = std::min(curlevel, maxlevelcopy); level >= 0; level--) {
                if (level > maxlevelcopy || level < 0)  // possible?
                    throw std::runtime_error("Level error");

                std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates = searchBaseLayer(
                        currObj, window_idx, level);

                currObj = mutuallyConnectNewElement(window_idx, top_candidates, level, false);
            }
        } else {
            // Do nothing for the first element
            enterpoint_node_ = window_idx;
            maxlevel_ = curlevel;
        }

        // Releasing lock for the maximum level
        if (curlevel > maxlevelcopy) {
            enterpoint_node_ = window_idx;
            maxlevel_ = curlevel;
        }
    }
};