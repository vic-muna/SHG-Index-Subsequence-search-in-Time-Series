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
        
        
        size_t link_list_size = sizeof(unsigned short) + M0 * sizeof(int);
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
        
        std::cout << "Level 0 memory allocated at: " << (void*) level0_memory_ << " | Node size: " << size_data_per_element_ << " bytes." << std::endl;
    }

    //TODO
    ~TimeSeriesHNSW(){
        if (level0_memory_ != nullptr) free(level0_memory_);
        if (linkLists_ != nullptr) {
            // Hier müssten später auch noch die individuellen Listen in linkLists_[i] gefree'd werden
            free(linkLists_);
        }
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
    
    void addDataPoint(size_t window_idx) {

        int curlevel = getRandomLevel(mult_);

        memset(level0_memory_ + window_idx * size_data_per_element_, 0, size_data_per_element_);
        memcpy(getDataByInternalId(window_idx), &window_idx, sizeof(size_t));

        std::vector<float> rep;
        rep.reserve(data_rep_size_ + window_size_);

        //TODO: fill in the get functions
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
      
        memcpy(data_rep_memory_+cur_c*data_rep_size_*sizeof(float), rep.data(), rep.size() * sizeof(float));

        per.timePreprocessing += t.elapsed();
        if(flag == -1){
            return cur_c;
        }
        //t.restart();
        if (curlevel) {
            linkLists_[cur_c] = (char *) malloc(size_links_per_element_ * curlevel + 1);
            if (linkLists_[cur_c] == nullptr)
                throw std::runtime_error("Not enough memory: addPoint failed to allocate linklist");
            memset(linkLists_[cur_c], 0, size_links_per_element_ * curlevel + 1);
        }

        if ((signed)currObj != -1) {
            if (curlevel < maxlevelcopy) {

                dist_t curdist = getDisByLevel(cur_c, currObj, 0);

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
                            dist_t d = getDisByLevel(cur_c, cand, 0);

                            if (d < curdist) {
                                curdist = d;
                                currObj = cand;
                                changed = true;
                            }
                        }
                    }
                }
            }

            bool epDeleted = isMarkedDeleted(enterpoint_copy);
            for (int level = std::min(curlevel, maxlevelcopy); level >= 0; level--) {
                if (level > maxlevelcopy || level < 0)  // possible?
                    throw std::runtime_error("Level error");

                std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates = searchBaseLayer(
                        currObj, data_point, label, level);
                if (epDeleted) {
                    //top_candidates.emplace(fstdistfunc_(data_point, getDataByInternalId(enterpoint_copy), dist_func_param_), enterpoint_copy);
                    top_candidates.emplace(getDisByLevel(label, enterpoint_copy, 0), enterpoint_copy);
                    if (top_candidates.size() > ef_construction_)
                        top_candidates.pop();
                }
                currObj = mutuallyConnectNewElement(data_point, cur_c, top_candidates, level, false);
            }
        } else {
            // Do nothing for the first element
            enterpoint_node_ = 0;
            maxlevel_ = curlevel;
        }

        // Releasing lock for the maximum level
        if (curlevel > maxlevelcopy) {
            enterpoint_node_ = cur_c;
            maxlevel_ = curlevel;
        }

        //per.timeBuildindex += t.elapsed();
        return cur_c;
    }
};