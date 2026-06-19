#include "../hnswlib/hnswlib.h"
#include <iostream>
#include <stdexcept>
#include <stdlib.h>
#include <memory>

template<typename dist_t>
class TimeSeriesHNSW : public hnswlib::AlgorithmInterface<dist_t>{
public:
    char *level0_memory_{nullptr};
    size_t size_data_per_element_{0};

    char **linkLists_{nullptr};

    TimeSeriesHNSW(
        hnswlib::SpaceInterface<dist_t> *space, 
        size_t window_size, 
        size_t num_windows, 
        size_t M = 16) {

        size_t M0 = M * 2;          // Anzahl der Kanten ist für das Basislevel (Level 0) am höchsten
        
        size_t link_list_size = sizeof(unsigned short) + M0 * sizeof(int);
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

        std::cout << "Level 0 memory allocated at: " << (void*) level0_memory_ << " | Node size: " << size_data_per_element_ << " bytes." << std::endl;
    }

    ~TimeSeriesHNSW(){
        if (level0_memory_ != nullptr) free(level0_memory_);
        if (linkLists_ != nullptr) {
            // Hier müssten später auch noch die individuellen Listen in linkLists_[i] gefree'd werden
            free(linkLists_);
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
    
};