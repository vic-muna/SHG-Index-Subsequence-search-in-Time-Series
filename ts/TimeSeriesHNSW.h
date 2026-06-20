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
    size_t data_rep_size_{0};
    unsigned int maxFixLevel_{0}; // base layer in HNSW is layer0
    unsigned int k_{4}; // for determining amount of layers

    char **linkLists_{nullptr};

    TimeSeriesHNSW(
        hnswlib::SpaceInterface<dist_t> *space, 
        size_t window_size, 
        size_t num_windows, 
        size_t M = 16) {

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

    //TODO
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
    
    void addDataPoint(const void *data_point, labeltype label, int flag, Performance &per, bool replace_deleted = false) {
        if ((allow_replace_deleted_ == false) && (replace_deleted == true)) {
            throw std::runtime_error("Replacement of deleted elements is disabled in constructor");
        }

        // lock all operations with element by label
        //std::unique_lock <std::mutex> lock_label(getLabelOpMutex(label));
        if (!replace_deleted) {
            addDataPoint(data_point, label, -1, flag, per);
            return;
        }
        // check if there is vacant place
        tableint internal_id_replaced;
        //std::unique_lock <std::mutex> lock_deleted_elements(deleted_elements_lock);
        bool is_vacant_place = !deleted_elements.empty();
        if (is_vacant_place) {
            internal_id_replaced = *deleted_elements.begin();
            deleted_elements.erase(internal_id_replaced);
        }
        //lock_deleted_elements.unlock();

        // if there is no vacant place then add or update point
        // else add point to vacant place
        if (!is_vacant_place) {
            addDataPoint(data_point, label, -1, flag, per);
        } else {
            // we assume that there are no concurrent operations on deleted element
            labeltype label_replaced = getExternalLabel(internal_id_replaced);
            setExternalLabel(internal_id_replaced, label);

            //std::unique_lock <std::mutex> lock_table(label_lookup_lock);
            label_lookup_.erase(label_replaced);
            label_lookup_[label] = internal_id_replaced;
            //lock_table.unlock();

            unmarkDeletedInternal(internal_id_replaced);
            updatePoint(data_point, internal_id_replaced, 1.0);
        }
    }
};