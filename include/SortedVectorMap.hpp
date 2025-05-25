#include <vector>
#include <algorithm>
#include <iostream>
#include <optional>

template<typename Key, typename Value>
class SortedVectorMap {
public:
    using Pair = std::pair<Key, Value>;

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

    // Find (const)
    const Value* find(const Key& key) const {
        auto it = lower_bound(key);
        if (it != data.end() && it->first == key) {
            return &it->second;
        }
        return nullptr;
    }

    // Find (non-const)
    Value* find(const Key& key) {
        auto it = lower_bound(key);
        if (it != data.end() && it->first == key) {
            return &it->second;
        }
        return nullptr;
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
};

