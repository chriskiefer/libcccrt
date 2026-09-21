// libcccrt core - ring buffer
// Dependency-free (no Eigen, no heap). Suitable for embedded targets.
#pragma once

#include <cstddef>

namespace cccrt {

// Fixed-capacity ring buffer over caller-owned storage.
// Push samples one at a time, then copy out the most recent `winSize` samples.
template <typename T>
class RingBuffer {
public:
    RingBuffer() = default;
    RingBuffer(T* storage, size_t capacity) { init(storage, capacity); }

    // `storage` must hold `capacity` elements and outlive the buffer. Zero-fills it.
    void init(T* storage, size_t capacity) {
        buf_ = storage;
        cap_ = capacity;
        idx_ = 0;
        for (size_t i = 0; i < cap_; ++i) buf_[i] = T(0);
    }

    void push(T x) {
        buf_[idx_] = x;
        if (++idx_ == cap_) idx_ = 0;
    }

    size_t capacity() const { return cap_; }

    // Copies the `winSize` most recent samples into `out` (oldest first), ending
    // `offset` samples before the most recently pushed one.
    // winSize + offset must not exceed capacity.
    void copyLatest(T* out, size_t winSize, size_t offset = 0) const {
        const size_t back = winSize + offset;
        size_t start = idx_ >= back ? idx_ - back : cap_ - (back - idx_);
        for (size_t i = 0; i < winSize; ++i) {
            out[i] = buf_[start];
            if (++start == cap_) start = 0;
        }
    }

private:
    T* buf_ = nullptr;
    size_t cap_ = 0;
    size_t idx_ = 0;
};

} // namespace cccrt
