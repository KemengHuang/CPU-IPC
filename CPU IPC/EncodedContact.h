#pragma once

#include <array>

// Compact legacy contact representation used as an ordered map key.
// Negative entries encode PP/PE/PT variants and multiplicity; decoding rules
// are centralized in the contact mechanics and broad-phase modules.
class EncodedContact {
public:
    explicit EncodedContact(int first = -1)
        : data_{ { first, -1, -1, -1 } }
    {
    }

    EncodedContact(int first, int second, int third, int fourth)
        : data_{ { first, second, third, fourth } }
    {
    }

    void set(int first, int second, int third, int fourth)
    {
        data_ = { { first, second, third, fourth } };
    }

    int& operator[](int index) { return data_[index]; }
    const int& operator[](int index) const { return data_[index]; }

    friend bool operator==(const EncodedContact& lhs, const EncodedContact& rhs)
    {
        return lhs.data_ == rhs.data_;
    }

    friend bool operator<(const EncodedContact& lhs, const EncodedContact& rhs)
    {
        return lhs.data_ < rhs.data_;
    }

private:
    std::array<int, 4> data_;
};
