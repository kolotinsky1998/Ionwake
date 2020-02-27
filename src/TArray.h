#ifndef OPTIMIZATION_TARRAY_H
#define OPTIMIZATION_TARRAY_H

#include <cstddef>
#include <array>
#include <vector>

template<class T, const std::size_t N>
class TArray {
public:

    TArray();

    explicit TArray(std::array<std::size_t, N>);

    TArray(TArray<double, N> &array);

    ~TArray();

    TArray &operator=(TArray<T, N> &&) noexcept;

    inline T &operator[](std::array<std::size_t, N>);

    inline T &operator[](size_t);

    void fill(const T &);

    std::size_t getTotalSize() noexcept;

    bool isEmpty() noexcept;

    bool swap(TArray<T, N> &other);

private:
    const std::array<std::size_t, N> sizes;
    T *data;
    std::array<std::size_t, N> shifts;
    std::size_t totalSize;
};


template<class T, std::size_t N>
inline T &TArray<T, N>::operator[](std::array<std::size_t, N> position) {
    std::size_t convertedPosition = 0;
    for (int i = 0; i < N; ++i) {
        convertedPosition += position[i] * shifts[i];
    }
    return data[convertedPosition];
}

template<class T, std::size_t N>
TArray<T, N>::TArray(std::array<std::size_t, N> sizes): sizes(sizes) {
    totalSize = 1;
    for (const auto &size : sizes) {
        totalSize *= size;
    }
    shifts[N - 1] = 1;
    for (int i = N - 2; i >= 0; --i) {
        shifts[i] = shifts[i + 1] * sizes[i + 1];
    }
    data = new T[totalSize];
}

template<class T, std::size_t N>
std::size_t TArray<T, N>::getTotalSize() noexcept {
    return totalSize;
}

template<class T, std::size_t N>
bool TArray<T, N>::isEmpty() noexcept {
    return N == 0 || totalSize == 0;
}

template<class T, std::size_t N>
void TArray<T, N>::fill(const T &value) {
    std::fill(data, data + totalSize, value);
}

template<class T, std::size_t N>
T &TArray<T, N>::operator[](size_t position) {
    return data[position];
}

template<class T, std::size_t N>
TArray<T, N>::~TArray() {
    delete[] data;
}

template<class T, std::size_t N>
bool TArray<T, N>::swap(TArray<T, N> &other) {
    if (this->totalSize != other.totalSize || this->sizes != other.sizes) {
        return false;
    }
    T *otherData = other.data;
    other.data = this->data;
    this->data = otherData;
    return true;
}

template<class T, std::size_t N>
TArray<T, N> &TArray<T, N>::operator=(TArray<T, N> &&other) noexcept {
    delete[] data;
    new(this) TArray(other);
    return *this;
}

template<class T, std::size_t N>
TArray<T, N>::TArray(TArray<double, N> &other):
        sizes(other.sizes),
        totalSize(other.totalSize),
        shifts(other.shifts) {
    data = new T[totalSize];
    std::copy(other.data, other.data + other.totalSize, data);
}

template<class T, std::size_t N>
TArray<T, N>::TArray() : sizes(std::array<std::size_t, N>()) {
    totalSize = 0;
    shifts = std::array<std::size_t, N>();
    data = nullptr;
}

#endif