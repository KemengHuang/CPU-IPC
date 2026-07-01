#pragma once

#include <cassert>
#include <cstdint>

class mesh3D;

class Triangle {
public:
    uint64_t key[3];

    Triangle(const uint64_t* p_key) {
        key[0] = p_key[0];
        key[1] = p_key[1];
        key[2] = p_key[2];
    }

    Triangle(uint64_t key0, uint64_t key1, uint64_t key2) {
        key[0] = key0;
        key[1] = key1;
        key[2] = key2;
    }

    uint64_t operator[](int i) const {
        assert(0 <= i && i <= 2);
        return key[i];
    }

    bool operator<(const Triangle& right) const {
        if (key[0] < right.key[0]) {
            return true;
        } else if (key[0] == right.key[0]) {
            if (key[1] < right.key[1]) {
                return true;
            } else if (key[1] == right.key[1]) {
                if (key[2] < right.key[2]) {
                    return true;
                }
            }
        }
        return false;
    }

    bool operator==(const Triangle& right) const {
        return key[0] == right[0] && key[1] == right[1] && key[2] == right[2];
    }
};

namespace MeshGeometry {

void computeSurface(mesh3D& mesh);

} // namespace MeshGeometry
