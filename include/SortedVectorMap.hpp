#pragma once
#include <vector>
#include <stdexcept>
#include <algorithm>


template<typename Key, typename Value>
class SortedVectorMap {
    using Pair = std::pair<Key, Value>;
private:
    std::vector<Pair> data;

    // Binary search lower bound
    typename std::vector<Pair>::iterator lower_bound(const Key& key) {
        return std::lower_bound(data.begin(), data.end(), key,
            [](const Pair& p, const Key& k) { return p.first < k; });
    }

    typename std::vector<Pair>::const_iterator lower_bound(const Key& key) const {
        return std::lower_bound(data.begin(), data.end(), key,
            [](const Pair& p, const Key& k) { return p.first < k; });
    }
public:

    // Insert or update
    void insert(const Key& key, const Value& value) {
        auto it = lower_bound(key);
        if (it != data.end() && it->first == key) {
            it->second = value; // update
        } else {
            data.insert(it, {key, value}); // insert in sorted order
        }
    }

    // Access or insert default value
    Value& operator[](const Key& key) {
        auto it = lower_bound(key);
        if (it != data.end() && it->first == key) {
            return it->second;
        } else {
            it = data.insert(it, {key, Value{}});
            return it->second;
        }
    }


    Value& at(const Key& key) {
        auto it = lower_bound(key);
#ifndef NDEBUG
		if (it == data.end()){
			throw std::out_of_range("Key not in structure");
		}
#endif
		return it->second; 
    }

    // Find (const)
	auto find(const Key& key) const {
        auto it = lower_bound(key);
        if (it != data.end() && it->first == key) {
            return it;
        }
        return data.end();
    }

    // Find (non-const)
    auto find(const Key& key) {
        auto it = lower_bound(key);
        if (it != data.end() && it->first == key) {
            return it;
        }
        return data.end();
    }

    // Erase by key
    bool erase(const Key& key) {
        auto it = lower_bound(key);
        if (it != data.end() && it->first == key) {
            data.erase(it);
            return true;
        }
        return false;
    }

	using iterator = typename std::vector<Pair>::iterator;
	using const_iterator = typename std::vector<Pair>::const_iterator;

    // Erase by iterator
    iterator erase(std::vector<Pair>::iterator it) {
		return data.erase(it);
    }

    // Size
    size_t size() const { return data.size(); }

    // Empty
    bool empty() const { return data.empty(); }

    // Iterators
    auto begin() { return data.begin(); }
    auto end() { return data.end(); }
    auto begin() const { return data.begin(); }
    auto end() const { return data.end(); }

    auto cbegin() const { return data.cbegin(); }
    auto cend() const { return data.cend(); }


};




/*
template<
std::convertible_to<std::size_t> Key,
	typename Value
	>
requires std::is_pointer_v<Value>
class FilteredVector {
private:
    std::vector<Value> data;
	std::size_t i_;

	void skip_nulls() {
		while (i_ < data.size() && data[i_] == nullptr) {
			++i_;
		}
	}

public:
    inline void insert(const Key& key, const Value& value) {
		data[key] = value;
    }

    // Access or insert default value
    inline Value& operator[](const Key& key) {
		return data[key];
    }


    inline Value& at(const Key& key) {
#ifndef NDEBUG
		if (data[key] == nullptr){
			throw std::out_of_range("Key not in structure");
		}
#endif
    }

    // Find (const)
    const Value* find(const Key& key) const {
		return &data[key];
    }

    // Find (non-const)
    Value* find(const Key& key) {
        return &data[key];
    }

    // Erase by key
    bool erase(const Key& key) {
		if (data[key] == nullptr){
			return false;
		}
		data[key]=nullptr;
		return true;
    }

    // Size
    size_t size() const { return data.size(); }

    // Empty
    bool empty() const { return data.empty(); }

    // Iterators
	//
	// Iterator
	template <typename VecPtr>
    class base_iterator {
    public:
        using value_type = std::pair<std::size_t, Value>;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::forward_iterator_tag;

        base_iterator(VecPtr vec, std::size_t i) : vec_(vec), i_(i) {
            skip_nulls();
        }

        value_type operator*() const {
            return {i_, (*vec_)[i_]};
        }

        base_iterator& operator++() {
            ++i_;
            skip_nulls();
            return *this;
        }

        bool operator==(const base_iterator& other) const {
            return i_ == other.i_;
        }

        bool operator!=(const base_iterator& other) const {
            return !(*this == other);
        }

    private:
        VecPtr vec_;
        std::size_t i_;

        void skip_nulls() {
            while (i_ < vec_->size() && (*vec_)[i_] == nullptr) {
                ++i_;
            }
        }
    };

	using iterator = base_iterator<std::vector<Value>*>;
	using const_iterator = base_iterator<const std::vector<Value>*>;
	 iterator begin() {
        return iterator(&data, 0);
    }

    iterator end() {
        return iterator(&data, data.size());
    }

    const_iterator begin() const {
        return const_iterator(&data, 0);
    }

    const_iterator end() const {
        return const_iterator(&data, data.size());
    }

    const_iterator cbegin() const {
        return begin();
    }

    const_iterator cend() const {
        return end();
    }

};
*/

template<
    std::convertible_to<std::size_t> Key,
    typename Value
>
requires std::is_pointer_v<Value>
class FilteredVector : public std::ranges::view_base
{
private:
    std::vector<Value> data;

public:
    inline void insert(const Key& key, const Value& value) {
        if (static_cast<std::size_t>(key) >= data.size()) {
            data.resize(static_cast<std::size_t>(key) + 1, nullptr);
        }
        data[key] = value;
    }

    inline Value& operator[](const Key& key) {
        if (static_cast<std::size_t>(key) >= data.size()) {
            data.resize(static_cast<std::size_t>(key) + 1, nullptr);
        }
        return data[key];
    }

    inline Value& at(const Key& key) {
#ifndef NDEBUG
        if (data[key] == nullptr){
            throw std::out_of_range("Key not in structure");
        }
#endif
        return data[key];
    }

    const Value* find(const Key& key) const {
        return &data[key];
    }

    Value* find(const Key& key) {
        return &data[key];
    }

    bool erase(const Key& key) {
        if (data[key] == nullptr) {
            return false;
        }
        data[key] = nullptr;
        return true;
    }

    std::size_t size() const { return data.size(); }

    bool empty() const { return data.empty(); }

    // Iterator
    template <typename VecPtr>
    class base_iterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = std::pair<std::size_t, Value>;
        using difference_type = std::ptrdiff_t;
        using pointer = void;
        using reference = value_type;

        base_iterator() = default;

        base_iterator(VecPtr vec, std::size_t i) : vec_(vec), i_(i) {
            skip_nulls();
        }

        value_type operator*() const {
            return {i_, (*vec_)[i_]};
        }

        base_iterator& operator++() {
            ++i_;
            skip_nulls();
            return *this;
        }

        bool operator==(const base_iterator& other) const {
            return i_ == other.i_;
        }

        bool operator!=(const base_iterator& other) const {
            return !(*this == other);
        }

    private:
        VecPtr vec_ = nullptr;
        std::size_t i_ = 0;

        void skip_nulls() {
            while (i_ < vec_->size() && (*vec_)[i_] == nullptr) {
                ++i_;
            }
        }
    };

    using iterator = base_iterator<std::vector<Value>*>;
    using const_iterator = base_iterator<const std::vector<Value>*>;

    iterator begin() {
        return iterator(&data, 0);
    }

    iterator end() {
        return iterator(&data, data.size());
    }

    const_iterator begin() const {
        return const_iterator(&data, 0);
    }

    const_iterator end() const {
        return const_iterator(&data, data.size());
    }

    const_iterator cbegin() const {
        return begin();
    }

    const_iterator cend() const {
        return end();
    }
};
