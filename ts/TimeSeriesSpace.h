#include "../hnswlib/hnswlib.h"
#include <vector>
#include <cmath>
#include <limits>

class TimeSeriesSpace : public hnswlib::SpaceInterface<float>{
private:
    float* data_;
    size_t window_size_;
    size_t num_windows_;

    std::vector<float> means_;
    std::vector<float> stds_;

    hnswlib::DISTFUNC<float> fstdistfunc_;

public:
    TimeSeriesSpace(const float* data, size_t window_size, size_t num_total_floats) {
        data_ = data;
        window_size_ = window_size;
        num_windows_ = num_total_floats - window_size + 1;

        means_.resize(num_windows_);
        stds_.resize(num_windows_);
        
        // TODO: calculate rolling means and standard deviations

        fstdistfunc_ = distance_wrapper;
    }   

    size_t get_data_size() override{
        return sizeof(size_t);
    }

    hnswlib::DISTFUNC<float> get_dist_func() override{
        return fstdistfunc_;
    }

    void *get_dist_func_param() override{
        return this;
    }

    ~TimeSeriesSpace() override{}


    static float distance_wrapper(const void* pVect1v, const void* pVect2v, const void* qty_ptr) {
        size_t idx1 = reinterpret_cast<const size_t*>(pVect1v);
        size_t idx2 = reinterpret_cast<const size_t*>(pVect2v);

        const TimeSeriesSpace* space = reinterpret_cast<const TimeSeriesSpace*>(qty_ptr);

        if (std::abs(static_cast<long long>(idx1)) - std::abs(static_cast<long long>(idx2)) <= space->window_size_ / 2){
            return std::numeric_limits<floats>::max();
        }

        const float* vec1 = space->data_ + idx1;
        const float* vec2 = space->data_ + idx2;
        float mu1 = space->means_[idx1];
        float sigma1 = space->stds_[idx1];
        float mu2 = space->means_[idx2];
        float sigma2 = space->stds_[idx2];

        float res = 0.0f;
        for (size_t i = 0; i < space->window_size_; ++i) {
            float norm1 = (vec1 - mu1) / sigma1;
            float norm2 = (vec2 - mu2) / sigma2;

            float diff = norm1 - norm2;
            res += diff * diff; 
        }

        return res;
    }
};