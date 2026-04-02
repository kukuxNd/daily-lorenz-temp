/*
 * Lorenz Attractor - 混沌理论的蝴蝶效应
 * 
 * 洛伦兹吸引子是最著名的混沌系统之一
 * 由气象学家 Edward Norton Lorenz 在 1963 年研究大气对流时发现
 * 
 * 方程组:
 *   dx/dt = σ(y - x)
 *   dy/dt = x(ρ - z) - y
 *   dz/dt = xy - βz
 * 
 * 典型参数: σ=10, ρ=28, β=8/3
 * 
 * 编译: g++ -O2 -o lorenz lorenz.cpp
 * 运行: ./lorenz
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>

// 洛伦兹方程参数
const double SIGMA = 10.0;
const double RHO   = 28.0;
const double BETA  = 8.0 / 3.0;

// 时间步长和模拟时长
const double DT = 0.01;
const int STEPS = 10000;

// 洛伦兹方程导数
void lorenz(double t, const double x[3], double dx[3]) {
    dx[0] = SIGMA * (x[1] - x[0]);
    dx[1] = x[0] * (RHO - x[2]) - x[1];
    dx[2] = x[0] * x[1] - BETA * x[2];
}

// 四阶龙格-库塔积分
void rk4(double& t, double x[3], double dt) {
    double k1[3], k2[3], k3[4], xm[3], xn[3];
    
    // k1
    lorenz(t, x, k1);
    for (int i = 0; i < 3; i++) xm[i] = x[i] + 0.5 * dt * k1[i];
    
    // k2
    lorenz(t + 0.5 * dt, xm, k2);
    for (int i = 0; i < 3; i++) xm[i] = x[i] + 0.5 * dt * k2[i];
    
    // k3
    lorenz(t + 0.5 * dt, xm, k3);
    for (int i = 0; i < 3; i++) xn[i] = x[i] + dt * k3[i];
    
    // k4
    double k4[3];
    lorenz(t + dt, xn, k4);
    
    // 合并
    for (int i = 0; i < 3; i++) {
        x[i] += (dt / 6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }
    t += dt;
}

// 投影到 XZ 平面做 ASCII 可视化
void asciiXZ(const std::vector<double>& xs, const std::vector<double>& zs) {
    const int W = 80, H = 40;
    char canvas[H][W];
    
    // 找到数据范围
    double xMin = 1e9, xMax = -1e9, zMin = 1e9, zMax = -1e9;
    for (size_t i = 0; i < xs.size(); i++) {
        if (xs[i] < xMin) xMin = xs[i];
        if (xs[i] > xMax) xMax = xs[i];
        if (zs[i] < zMin) zMin = zs[i];
        if (zs[i] > zMax) zMax = zs[i];
    }
    
    // 初始化空白画布
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
            canvas[y][x] = ' ';
    
    // 绘制轨迹
    for (size_t i = 1; i < xs.size(); i++) {
        int px = int((xs[i] - xMin) / (xMax - xMin) * (W - 1));
        int py = int((zs[i] - zMin) / (zMax - zMin) * (H - 1));
        py = H - 1 - py; // 反转Y轴
        if (px >= 0 && px < W && py >= 0 && py < H) {
            if (canvas[py][px] == ' ') {
                // 密度着色
                int density = int(8.0 * i / xs.size());
                canvas[py][px] = " .:;+*#@" [std::min(density, 7)];
            }
        }
    }
    
    // 输出
    std::cout << "\n╭" << std::string(W, '─') << "╮\n";
    for (int y = 0; y < H; y++) {
        std::cout << "│";
        for (int x = 0; x < W; x++)
            std::cout << canvas[y][x];
        std::cout << "│\n";
    }
    std::cout << "╰" << std::string(W, '─') << "╯\n";
    std::cout << "  X ∈ [" << std::fixed << std::setprecision(1) << xMin << ", " << xMax << "]\n";
    std::cout << "  Z ∈ [" << zMin << ", " << zMax << "]\n";
}

int main() {
    std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
    std::cout << "   🦋 Lorenz Attractor 洛伦兹吸引子\n";
    std::cout << "   混沌理论的蝴蝶效应 · 最小核心实现\n";
    std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
    std::cout << "\n参数: σ=" << SIGMA << " ρ=" << RHO << " β=" << BETA << "\n";
    std::cout << "步长: " << DT << " | 迭代: " << STEPS << " 步\n\n";
    
    // 两组初始条件，仅差 0.0001 —— 演示蝴蝶效应
    std::vector<double> x1(3), x2(3);
    x1[0] = 1.0; x1[1] = 1.0; x1[2] = 1.0;
    x2[0] = 1.0001; x2[1] = 1.0; x2[2] = 1.0;  // 仅此处不同！
    
    std::vector<double> t1(STEPS), t2(STEPS);
    std::vector<double> xs1(STEPS), xs2(STEPS), zs1(STEPS), zs2(STEPS);
    
    double t = 0.0;
    
    std::cout << "正在计算轨迹 (红色=初始[1,1,1], 蓝色=初始[1.0001,1,1])...\n\n";
    
    for (int i = 0; i < STEPS; i++) {
        rk4(t, x1.data(), DT);
        xs1[i] = x1[0]; zs1[i] = x1[2];
    }
    
    t = 0.0;
    for (int i = 0; i < STEPS; i++) {
        rk4(t, x2.data(), DT);
        xs2[i] = x2[0]; zs2[i] = x2[2];
    }
    
    // ASCII 可视化
    std::cout << "XZ 平面投影 (蝴蝶形状):\n";
    asciiXZ(xs1, zs1);
    
    // 蝴蝶效应量化
    double divMax = 0.0;
    for (int i = 0; i < STEPS; i++) {
        double d = std::sqrt(
            (xs1[i]-xs2[i])*(xs1[i]-xs2[i]) +
            (zs1[i]-zs2[i])*(zs1[i]-zs2[i])
        );
        if (d > divMax) divMax = d;
    }
    
    std::cout << "\n🦋 蝴蝶效应分析:\n";
    std::cout << "   初始距离: 0.0001\n";
    std::cout << "   最大偏离: " << std::scientific << divMax << "\n";
    std::cout << "   放大倍数: " << std::fixed << divMax / 0.0001 << "x\n\n";
    
    std::cout << "混沌特征: 对初值敏感 —— 毫厘之差，终成天壤之别\n";
    
    // 输出 gnuplot 数据文件
    std::ofstream f("lorenz.dat");
    for (int i = 0; i < STEPS; i++)
        f << xs1[i] << " " << xs1[i] << " " << zs1[i] << "\n";
    f.close();
    std::cout << "\n已生成 lorenz.dat (可用于 gnuplot 3D 绘图)\n";
    
    std::cout << "\n# gnuplot 3D 绘图命令:\n";
    std::cout << "# splot 'lorenz.dat' with lines notitle\n";
    std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
    
    return 0;
}
