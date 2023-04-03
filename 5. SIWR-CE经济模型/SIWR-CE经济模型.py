import math
import numpy as np
import matplotlib.pyplot as plt

N = 100000              # 总人口数
i0 = 1                  # 输入感染人数
r0 = 2                 # R0值，即一个传染源一天能传染多少人
r = 15                  # 恢复天数
Gc = 1000               # 单人单日贡献值
Ge = 1000               # 单人单日消耗值
s = 9                   # 一个大白管几个人
Day = 100               # 观察期长度

I = i0                 # I: illness     -> 得病者，初始值为输入感染人数
W = math.ceil(I / s)   # W: white       -> 大白
R = 0                  # R: recover     -> 恢复者
S = N - I - W - R      # S: susceptible -> 易感者，即易感人群

economy_list = []

I_list = [I]
W_list = [W]
R_list = [R]
S_list = [S]

ill_dict = dict(zip(range(1, r + 1, 1), [0] * 15))
ill_dict[1] = 1

def ill_dis():
    global ill_dict
    recover = ill_dict[r]
    for i in range(r - 1, 0, -1):
        ill_dict[i + 1] = ill_dict[i]
    return recover

def cal_per_day(day):
    global I, W, R, S, ill_dict, I_list, W_list, R_list, E_list
    It = I * r0

    max_e = math.ceil(s / (s + 1) * S)        # 最大的可能新增感染数
    It = max_e if It > max_e else It
    Rt = ill_dis()

    I = I + It - Rt
    W = math.ceil(I / s)
    R = R + Rt
    S = N - I - W - R
    ill_dict[1] = It

    I_list.append(I)
    W_list.append(W)
    R_list.append(R)
    S_list.append(S)

    print('day', day, I, W, R, S, I+W+R+S, max_e, ill_dict)
    return I, W, R, S

def cal_economy(I, W, R, S):
    global alleconomy
    creation = (S + R) * Gc
    expenditrue = (I + W) * Ge
    economy_list.append(creation - expenditrue)

def draw_EIRW():
    t = np.arange(0, Day + 1, 1)
    plt.figure(figsize=(16, 9))
    plt.plot(t, S_list, '-', label='S_list')
    plt.plot(t, I_list, '-', label='I_list')
    plt.plot(t, W_list, '-',  label='W_list')
    plt.plot(t, R_list, '-',  label='R_list')
    plt.legend()
    # plt.show()
    plt.savefig('./r0为20的模型效果.jpg')

def draw_economy():
    t = np.arange(0, Day, 1)
    plt.figure(figsize=(16, 9))
    plt.plot(t, economy_list, '-')
    # plt.show()
    plt.savefig('./r0为2的经济效果.jpg')

if __name__ == '__main__':
    for i in range(1, Day + 1):
        I, W, R, E = cal_per_day(i)
        cal_economy(I, W, R, S)
    # draw_EIRW()
    draw_economy()
    print(sum(economy_list))