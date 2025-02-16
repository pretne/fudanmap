#include <iostream>
#include <sstream>
#include <unordered_set>
#include <map>
#include <vector>
#include "tinyxml/tinyxml.h"  // 用于加载和解析OpenStreetMap的XML数据
#include "jsoncpp/json/json.h" // 用于JSON数据格式化和处理
#include "algorithm/Graph.h"   // 导入Graph类，包含图的表示和最短路径计算的实现
#include <chrono>               // 用于时间测量
#include <emscripten/bind.h>   // emscripten的绑定库，用于将C++函数暴露给JavaScript调用

using namespace emscripten;
using namespace std;

// 用来存储节点和ways的JSON数据
Json::Value nodes;
Json::Value ways;

// 创建一个图对象，graph用于存储所有的节点、边等图数据
Graph graph;

// 加载并解析XML文件，构建图
void load()
{
    // 加载地图文件
    TiXmlDocument tinyXmlDoc("map");
    tinyXmlDoc.LoadFile();
    TiXmlElement *root = tinyXmlDoc.RootElement();

    // 解析所有的节点元素
    TiXmlElement *nodeElement = root->FirstChildElement("node");
    for (; nodeElement; nodeElement = nodeElement->NextSiblingElement("node"))
    {
        Json::Value node;
        node["id"] = nodeElement->Attribute("id");  // 获取节点ID
        node["lon"] = nodeElement->Attribute("lon");  // 获取经度
        node["lat"] = nodeElement->Attribute("lat");  // 获取纬度

        // 将节点添加到图中
        graph.addNode(node["id"].asString(), stod(node["lat"].asString()), stod(node["lon"].asString()));

        // 保存节点数据到JSON对象
        nodes[nodeElement->Attribute("id")] = node;
    }

    // 解析所有的way元素
    TiXmlElement *wayElement = root->FirstChildElement("way");
    for (; wayElement; wayElement = wayElement->NextSiblingElement("way"))
    {
        Json::Value way;
        way["id"] = wayElement->Attribute("id");  // 获取way的ID

        Json::Value wayNodes;  // 用于存储与way相关的所有节点
        TiXmlElement *childNode = wayElement->FirstChildElement("nd");
        string firstNode = childNode->Attribute("ref");  // 获取第一个节点ID
        string lastNode = firstNode;  // 默认第一个节点和最后一个节点相同
        for (; childNode; childNode = childNode->NextSiblingElement("nd"))
        {
            string ref = childNode->Attribute("ref");  // 获取每个节点ID
            lastNode = ref;
            wayNodes.append(nodes[ref]);  // 将节点ID添加到wayNodes中
        }

        // 如果第一个节点和最后一个节点相同，认为这是一个环，跳过
        if (firstNode == lastNode)
        {
            continue;
        }

        way["nodes"] = wayNodes;

        // 使用KD树添加连接的边
        for (int i = 0; i < wayNodes.size() - 1; i++)
        {
            string id1 = wayNodes[i]["id"].asString();
            string id2 = wayNodes[i + 1]["id"].asString();
            graph.addEdge(id1, id2);  // 为每对节点添加边
        }

        Json::Value wayTags;
        TiXmlElement *childTag = wayElement->FirstChildElement("tag");
        for (; childTag; childTag = childTag->NextSiblingElement("tag"))
        {
            string name = childTag->Attribute("k");  // 获取标签的键
            string value = childTag->Attribute("v");  // 获取标签的值
            wayTags[name] = value;
        }
        way["tags"] = wayTags;

        // 将way数据添加到JSON对象中
        ways[wayElement->Attribute("id")] = way;
    }

    // 存储待删除的孤立节点
    std::vector<std::string> nodesToRemove;

    // 遍历所有节点，找出孤立的节点
    for (auto &node : nodes)
    {
        // 如果该节点没有邻居节点，说明是孤立节点
        if (graph.getNeighborsSize(node["id"].asString()) == 0)
        {
            nodesToRemove.push_back(node["id"].asString());  // 记录孤立节点ID
        }
    }

    // 删除孤立节点和相关的JSON数据
    for (const auto &nodeId : nodesToRemove)
    {
        graph.removeNode(nodeId);  // 从图中移除该孤立节点
        nodes.removeMember(nodeId);  // 从JSON数据中移除该孤立节点
    }
}

// 获取所有节点的JSON数据
string getNodes()
{
    Json::StreamWriterBuilder builder;
    string s = Json::writeString(builder, nodes);  // 将JSON对象转换为字符串
    return s;
}

// 获取所有ways的JSON数据
string getWays()
{
    Json::StreamWriterBuilder builder;
    string s = Json::writeString(builder, ways);  // 将JSON对象转换为字符串
    return s;
}

// 使用Dijkstra算法计算最短路径
string getShortestPathByDijkstra(string start, string end, string interruptedPointsJson)
{
    // 解析中断点信息
    Json::Value interruptedPointsJsonValue;
    Json::CharReaderBuilder readerBuilder;
    std::string errs;
    std::istringstream iss(interruptedPointsJson);  // 将传入的JSON字符串转换为流
    if (!Json::parseFromStream(readerBuilder, iss, &interruptedPointsJsonValue, &errs))
    {
        throw std::runtime_error("Failed to parse interruptedPointsJson: " + errs);  // 如果解析失败，抛出异常
    }

    // 将中断点保存为集合
    std::unordered_set<string> interruptedPoints;
    for (const auto &point : interruptedPointsJsonValue)
    {
        interruptedPoints.insert(point.asString());
    }

    // 记录开始时间
    auto startTime = std::chrono::high_resolution_clock::now();

    // 调用图的dijkstra算法，计算最短路径
    auto result = graph.dijkstra(start, end, interruptedPoints);

    // 记录结束时间
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = endTime - startTime;  // 计算算法执行时间

    // 构造返回的JSON对象，包含路径、距离和执行时间
    Json::Value path;
    for (const string &node : result.first)
    {
        path.append(nodes[node]);  // 将路径上的每个节点添加到JSON路径中
    }

    Json::Value res;
    res["path"] = path;
    res["distance"] = result.second;
    res["duration"] = duration.count();  // 记录算法执行时间

    Json::StreamWriterBuilder builder;
    string s = Json::writeString(builder, res);  // 将结果转为JSON字符串
    return s;
}

// 使用A*算法计算最短路径
string getShortestPathByAstar(string start, string end, string interruptedPointsJson)
{
    // 解析中断点信息
    Json::Value interruptedPointsJsonValue;
    Json::CharReaderBuilder readerBuilder;
    std::string errs;
    std::istringstream iss(interruptedPointsJson);  // 将传入的JSON字符串转换为流
    if (!Json::parseFromStream(readerBuilder, iss, &interruptedPointsJsonValue, &errs))
    {
        throw std::runtime_error("Failed to parse interruptedPointsJson: " + errs);  // 如果解析失败，抛出异常
    }

    // 将中断点保存为集合
    std::unordered_set<string> interruptedPoints;
    for (const auto &point : interruptedPointsJsonValue)
    {
        interruptedPoints.insert(point.asString());
    }

    // 记录开始时间
    auto startTime = std::chrono::high_resolution_clock::now();

    // 调用图的aStar算法，计算最短路径
    auto result = graph.aStar(start, end, interruptedPoints);

    // 记录结束时间
    auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = endTime - startTime;  // 计算算法执行时间

    // 构造返回的JSON对象，包含路径、距离和执行时间
    Json::Value path;
    for (const string &node : result.first)
    {
        path.append(nodes[node]);  // 将路径上的每个节点添加到JSON路径中
    }

    Json::Value res;
    res["path"] = path;
    res["distance"] = result.second;
    res["duration"] = duration.count();  // 记录算法执行时间

    Json::StreamWriterBuilder builder;
    string s = Json::writeString(builder, res);  // 将结果转为JSON字符串
    return s;
}

// emscripten绑定部分，暴露C++函数给JavaScript调用
EMSCRIPTEN_BINDINGS()
{
    emscripten::function("load", &load); // 绑定load函数
    emscripten::function("getNodes", &getNodes); // 绑定getNodes函数
    emscripten::function("getWays", &getWays); // 绑定getWays函数
    emscripten::function("getShortestPathByDijkstra", &getShortestPathByDijkstra); // 绑定getShortestPathByDijkstra函数
    emscripten::function("getShortestPathByAstar", &getShortestPathByAstar); // 绑定getShortestPathByAstar函数
}
