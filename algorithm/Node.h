#ifndef NODE_H
#define NODE_H

#include <string>
#include <cmath>

class Node {
private:
    std::string id;
    double lat; // 纬度
    double lon; // 经度

public:
    // 构造函数
    Node(const std::string& id, double lat, double lon) : id(id), lat(lat), lon(lon) {}

    // 获取ID
    std::string getId() const {
        return id;
    }

    // 获取纬度
    double getLat() const {
        return lat;
    }

    // 获取经度
    double getLon() const {
        return lon;
    }

    // 计算与另一节点的欧几里得距离
    double distance(const Node* other) const {
        double dLat = lat - other->getLat();
        double dLon = lon - other->getLon();
        return std::sqrt(dLat * dLat + dLon * dLon);
    }
};

#endif
