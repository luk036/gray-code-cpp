# Gray Code C++ 项目

## 项目概述

这是一个 C++ 实现的格雷码（Gray Code）算法项目，包含两个主要模块：

1. **middle** - 实现中间层格雷码算法，基于 [Muetze, Nummenpalo] 的研究
2. **rect** - 实现各种类型的矩形剖分（rectangulation）生成算法，基于 [Merino, Muetze] 的研究

项目使用 C++11 标准编写，采用 GNU GPL v2/v3 许可证。

## 项目结构

```
gray-code-cpp/
├── middle/                 # 中间层格雷码模块
│   ├── hamcycle.cpp/hpp    # Hamilton 循环实现
│   ├── tree.cpp/hpp        # 树结构实现
│   ├── vertex.cpp/hpp      # 顶点表示
│   ├── main.cpp            # 主程序入口
│   └── Makefile            # 构建配置
├── rect/                   # 矩形剖分模块
│   ├── rectangulation.cpp/hpp  # 矩形剖分核心实现
│   ├── edge.cpp/hpp        # 边表示
│   ├── rectangle.cpp/hpp   # 矩形表示
│   ├── vertex.cpp/hpp      # 顶点表示
│   ├── wall.cpp/hpp        # 墙表示
│   ├── main.cpp            # 主程序入口
│   └── Makefile            # 构建配置
├── .github/workflows/      # GitHub Actions 工作流
│   └── jekyll-gh-pages.yml # Jekyll 部署配置
└── .xmake/                 # XMake 构建系统缓存
```

## 构建和运行

### middle 模块

```bash
cd middle
make                    # 编译生成 middle 可执行文件
./middle -h             # 显示帮助信息
./middle -n2            # 计算长度为 5 的中间层格雷码
./middle -n2 -v01010    # 从指定位串开始计算
./middle -n2 -p1        # 打印翻转位置而非位串
./middle -n10 -l50      # 限制输出数量
```

**middle 参数说明：**
- `-n{1,2,...}`: 指定参数 n，生成长度为 2n+1、权重为 n 或 n+1 的位串
- `-l{-1,0,1,2,...}`: 要列出的位串数量，-1 表示完整循环
- `-v{0,1}^{2n+1}`: 初始位串（长度 2n+1，权重为 n 或 n+1）
- `-s{0,1}`: 是否存储并打印所有访问的位串
- `-p{0,1}`: 打印翻转位置而非位串

### rect 模块

```bash
cd rect
make                    # 编译生成 rect 可执行文件
./rect -h               # 显示帮助信息
./rect -n5 -c           # 生成 5 个矩形的剖分并计数
./rect -n5 -t2 -c       # 生成对角类型的矩形剖分
./rect -n5 -p3456 -c    # 生成避开特定模式的矩形剖分
./rect -n10 -t3 -l30    # 生成块对齐类型的矩形剖分，限制输出
```

**rect 参数说明：**
- `-n{1,2,...}`: 矩形数量
- `-t{1,2,3}`: 矩形剖分基础类型（1=通用，2=对角，3=块对齐）
- `-p{1,2,..,8}`: 禁止模式（1=顺时针风车，2=逆时针风车，3-8=其他模式）
- `-l{-1,0,1,2,...}`: 要列出的矩形剖分数量，-1 表示全部
- `-q`: 静默输出
- `-c`: 输出矩形数量

## 开发约定

### 编码标准
- 使用 C++11 标准
- 遵循 GNU GPL 许可证
- 代码包含详细的版权声明和许可证信息
- 使用 Doxygen 风格的注释

### 构建系统
- 使用传统 Makefile 构建系统
- 编译器：g++
- 编译标志：`-std=c++11 -O3`（middle 模块使用 C++0x）
- 静态链接：`-static`

### 测试和验证
- 每个模块都有独立的 main.cpp 作为测试入口
- 提供详细的命令行帮助信息
- 包含多种参数组合的示例用法

## 项目特点

1. **算法实现**: 实现了组合数学中的经典格雷码算法和矩形剖分生成算法
2. **模块化设计**: 两个独立的功能模块，各自有完整的实现
3. **性能优化**: 使用 O3 优化级别，可选择 NVISIT 宏进一步提升性能
4. **学术背景**: 基于最新的组合数学研究成果实现

## 许可证

本项目采用 GNU General Public License (GPL)：
- middle 模块：GPL v2
- rect 模块：GPL v3

## 贡献

项目由 Torsten Muetze、Jerri Nummenpalo 和 Arturo Merino 等人开发和维护。