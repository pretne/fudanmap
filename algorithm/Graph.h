#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <map>
#include <unordered_set>
#include <iostream>
#include <string>
#include <queue>
#include <limits>
#include <algorithm>
#include <cmath>
#include "Node.h" // 引入Node.h，Node类定义了节点的基本信息和计算距离的功能

using namespace std;

// 定义KD树节点
// KD树节点用于存储空间中的节点信息，并支持高效的最近邻搜索
struct KDNode {
    Node* node;   // 每个KD树节点包含一个Node指针，表示该节点
    KDNode* left; // 左子树
    KDNode* right; // 右子树

    // 构造函数，初始化节点和左右子树指针
    KDNode(Node* n) : node(n), left(nullptr), right(nullptr) {}
};

class Graph {
private:
    // 存储所有的节点，使用string作为ID
    map<string, Node*> nodes;

    // 邻接表，存储节点间的边及其权重
    map<string, map<string, double>> adjList;

    // KD树的根节点
    KDNode* kdRoot;

    // 构建KD树的辅助函数
    // 参数：points为节点列表，depth为当前递归的深度
    KDNode* buildKDTree(vector<Node*>& points, int depth) {
        if (points.empty()) return nullptr;  // 如果没有节点，返回空指针

        // 按照当前深度选择分割维度：纬度（depth % 2 == 0）或经度（depth % 2 == 1）
        int axis = depth % 2;
        // 按照当前维度对节点进行排序
        sort(points.begin(), points.end(), [axis](Node* a, Node* b) {
            return axis == 0 ? a->getLat() < b->getLat() : a->getLon() < b->getLon();
        });

        // 选择中位数节点作为当前节点
        int medianIndex = points.size() / 2;
        Node* medianNode = points[medianIndex];

        // 创建当前节点，并递归构建左右子树
        KDNode* node = new KDNode(medianNode);
        vector<Node*> leftPoints(points.begin(), points.begin() + medianIndex);  // 左半部分
        vector<Node*> rightPoints(points.begin() + medianIndex + 1, points.end());  // 右半部分

        node->left = buildKDTree(leftPoints, depth + 1);
        node->right = buildKDTree(rightPoints, depth + 1);

        return node;  // 返回当前树的根节点
    }

    // 最近邻搜索的辅助函数
    // 目标是找到距离目标节点最近的节点
    void nearestNeighborSearch(KDNode* root, const Node* target, int depth, 
                                KDNode*& best, double& bestDist) const {
        if (!root) return;  // 如果当前节点为空，结束递归

        // 计算当前节点与目标节点的距离
        double d = root->node->distance(target);
        if (d < bestDist) {
            bestDist = d;  // 更新最短距离
            best = root;   // 更新最近邻节点
        }

        // 计算当前深度对应的分割维度：纬度或经度
        int axis = depth % 2;
        // 计算目标节点与当前节点在分割维度上的差值
        double diff = (axis == 0) ? target->getLat() - root->node->getLat() 
                                  : target->getLon() - root->node->getLon();

        // 选择搜索左子树或右子树
        KDNode* near = (diff < 0) ? root->left : root->right;
        KDNode* far = (diff < 0) ? root->right : root->left;

        // 递归搜索较可能包含最近邻的子树
        nearestNeighborSearch(near, target, depth + 1, best, bestDist);

        // 如果分割面距离最近邻距离小于当前最佳距离，递归搜索另一子树
        if (abs(diff) < bestDist) {
            nearestNeighborSearch(far, target, depth + 1, best, bestDist);
        }
    }

public:
    // 构造函数，初始化kdRoot为nullptr
    Graph() : kdRoot(nullptr) {}

    // 析构函数，释放所有节点和KD树的内存
    ~Graph() {
        for (auto& node : nodes) {
            delete node.second;  // 删除每个节点
        }
        delete kdRoot;  // 删除KD树
    }

    // 添加一个节点，指定节点ID、纬度和经度
    void addNode(string id, double lat, double lon) {
        if (nodes.find(id) == nodes.end()) {
            Node* node = new Node(id, lat, lon);
            nodes[id] = node;  // 将新节点添加到图中
        } else {
            cout << "Node with id " << id << " already exists!" << endl;  // 节点已存在时提示
        }
    }

    // 构建KD树
    void buildTree() {
        vector<Node*> points;
        for (auto& pair : nodes) {
            points.push_back(pair.second);  // 将所有节点放入一个vector中
        }
        kdRoot = buildKDTree(points, 0);  // 构建KD树
    }

    // 查找距离指定位置最近的节点
    Node* findNearestNeighbor(double lat, double lon) const {
        if (!kdRoot) return nullptr;  // 如果KD树为空，返回nullptr

        Node target("", lat, lon);  // 创建目标节点
        KDNode* best = nullptr;
        double bestDist = numeric_limits<double>::infinity();  // 初始最佳距离为无穷大

        nearestNeighborSearch(kdRoot, &target, 0, best, bestDist);  // 执行最近邻搜索

        return best ? best->node : nullptr;  // 返回找到的最近邻节点
    }

    // 打印图中的所有节点信息
    void printGraph() {
        cout << "Graph nodes:" << endl;
        for (auto& node : nodes) {
            cout << "Node ID: " << node.first << " Location: (" 
                 << node.second->getLat() << ", " 
                 << node.second->getLon() << ")" << endl;
        }
    }

    // 获取节点数量
    int getNodeCount() const {
        return nodes.size();
    }

    // 计算两节点之间的距离
    double getDistance(string id1, string id2) {
        if (nodes.find(id1) != nodes.end() && nodes.find(id2) != nodes.end()) {
            return nodes[id1]->distance(nodes[id2]);  // 使用Node类的distance函数计算距离
        }
        return -1;  // 如果节点不存在，返回-1表示错误
    }

    // 删除一个节点
    void removeNode(string id) {
        if (nodes.find(id) != nodes.end()) {
            nodes.erase(id);  // 删除指定ID的节点
        } else {
            cout << "Node with id " << id << " does not exist!" << endl;  // 如果节点不存在，输出错误信息
        }
    }

    // 添加一条边，指定两节点ID
    void addEdge(string id1, string id2) {
        if (nodes.find(id1) != nodes.end() && nodes.find(id2) != nodes.end()) {
            double distance = nodes[id1]->distance(nodes[id2]);  // 计算两节点间的距离
            adjList[id1][id2] = distance;  // 添加到邻接表中
            adjList[id2][id1] = distance;  // 因为是无向图，所以反向也要添加
        } else {
            cout << "One or both nodes do not exist!" << endl;  // 如果节点不存在，输出错误信息
        }
    }

    // 获取指定节点的邻居数量
    int getNeighborsSize(string id) {
        if (adjList.find(id) != adjList.end()) {
            return adjList[id].size();
        }
        return 0;
    }

    // 使用Dijkstra算法计算最短路径
    pair<vector<string>, double> dijkstra(const string& start, const string& end, const std::unordered_set<std::string>& interruptedPoints) {
        // 初始化距离和前驱节点
        map<string, double> dist;
        map<string, string> prev;
        priority_queue<pair<double, string>, vector<pair<double, string>>, greater<>> pq;  // 最小堆

        // 初始化所有节点的距离为无穷大
        for (const auto& node : nodes) {
            dist[node.first] = numeric_limits<double>::infinity();
            prev[node.first] = "";
        }
        dist[start] = 0;  // 起始节点距离为0
        pq.push({0, start});  // 将起始节点加入优先队列

        while (!pq.empty()) {
            string u = pq.top().second;
            pq.pop();

            // 如果当前节点是中断点，跳过
            if (interruptedPoints.count(u)) continue;

            // 如果已到达终点，停止搜索
            if (u == end) break;

            // 遍历当前节点的所有邻居
            for (const auto& neighbor : adjList[u]) {
                string v = neighbor.first;

                // 如果邻居是中断点，跳过
                if (interruptedPoints.count(v)) continue;

                double weight = neighbor.second;
                double alt = dist[u] + weight;

                if (alt < dist[v]) {
                    dist[v] = alt;  // 更新最短距离
                    prev[v] = u;    // 更新前驱节点
                    pq.push({dist[v], v});  // 将邻居加入优先队列
                }
            }
        }

        // 构建最短路径
        vector<string> path;
        string curr = end;

        while (!curr.empty()) {
            path.push_back(curr);
            curr = prev[curr];  // 通过前驱节点回溯路径
        }

        // 如果路径为空，表示没有找到从起点到终点的路径
        if (path.size() == 1 && path[0] != start) return {{}, -1};

        reverse(path.begin(), path.end());  // 反转路径，得到从起点到终点的正确顺序
        return {path, dist[end]};  // 返回路径和最短距离
    }

    // 使用A*算法计算最短路径
    pair<vector<string>, double> aStar(const string& start, const string& end, const std::unordered_set<std::string>& interruptedPoints) {
        map<string, double> g_score;
        map<string, double> f_score;
        map<string, string> prev;
        priority_queue<pair<double, string>, vector<pair<double, string>>, greater<>> pq;

        // 初始化g_score和f_score
        for (auto& pair : nodes) {
            g_score[pair.first] = numeric_limits<double>::infinity();
            f_score[pair.first] = numeric_limits<double>::infinity();
        }

        // 起点初始化
        g_score[start] = 0;
        f_score[start] = heuristic(start, end);  // 计算启发式估价
        pq.push({f_score[start], start});  // 将起始节点加入队列

        while (!pq.empty()) {
            string current = pq.top().second;
            pq.pop();

            // 如果当前节点是中断点，跳过
            if (interruptedPoints.count(current)) continue;

            // 如果到达终点，停止搜索
            if (current == end) break;

            // 遍历当前节点的所有邻居
            for (auto& neighbor : adjList[current]) {
                string neighbor_id = neighbor.first;

                // 如果邻居是中断点，跳过
                if (interruptedPoints.count(neighbor_id)) continue;

                double weight = neighbor.second;
                double tentative_g_score = g_score[current] + weight;

                if (tentative_g_score < g_score[neighbor_id]) {
                    prev[neighbor_id] = current;
                    g_score[neighbor_id] = tentative_g_score;
                    f_score[neighbor_id] = g_score[neighbor_id] + heuristic(neighbor_id, end);
                    pq.push({f_score[neighbor_id], neighbor_id});
                }
            }
        }

        // 如果没有找到路径，返回空路径和-1
        if (prev.find(end) == prev.end()) return {{}, -1};

        // 构建路径
        vector<string> path;
        for (string at = end; !at.empty(); at = prev[at]) {
            path.push_back(at);
        }
        reverse(path.begin(), path.end());

        return {path, g_score[end]};  // 返回路径和最短路径的代价
    }

    // 计算启发式距离（从当前节点到目标节点的估算距离）
    double heuristic(const string& current, const string& end) {
        if (nodes.find(current) != nodes.end() && nodes.find(end) != nodes.end()) {
            Node* currentNode = nodes[current];
            Node* endNode = nodes[end];
            return currentNode->distance(endNode);  // 计算当前节点到目标节点的距离作为启发式估算
        }
        return 0.0;
    }
};

#endif
