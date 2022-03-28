'''
@coding: utf-8
@author: moyuweiqing
'''

import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def dySEIDR_C(y, t, alpha, beta, gamma, delta, lamda):
    s, e, i, d, r = y
    n = s + e + i + d + r

    if t >= 30 and t <= 60:                 # 从第30天开始病毒发生变异
        beta = beta * 1.2                   # 传染性增强
        delta = delta * 0.8                 # 致死性降低
    elif t > 60:                            # 第二次变异
        beta = beta * 1.5                   # 传染性进一步增强
        delta = delta * 0.5                 # 致死性进一步降低

    ds_dt = - lamda * beta * s * (i + e) / n
    de_dt = lamda * beta * s * (i + e) / n - alpha * e
    di_dt = alpha * e - gamma * i - delta * i
    dd_dt = delta * i
    dr_dt = gamma * i

    return np.array([ds_dt, de_dt, di_dt, dd_dt, dr_dt])


def draw_SEIDR(ySEIR):
    plt.figure(figsize=(16, 9))
    plt.plot(t, ySEIR[:, 0], '--', label='s(t)--SEIR')
    plt.plot(t, ySEIR[:, 1], '-.', label='e(t)--SEIR')
    plt.plot(t, ySEIR[:, 2], '-', label='i(t)--SEIR')
    plt.plot(t, ySEIR[:, 3], '.', label='d(t)--SEIR')
    plt.plot(t, ySEIR[:, 4], '--', label='r(t)--SEIR')
    plt.legend()
    # plt.show()
    plt.savefig('./img/' + '2.4 SEIDR_C_model.jpg')


if __name__ == '__main__':
    number = 10000  # 模型总人数

    s = 9990  # Susceptible 易感人群数量
    e = 0  # Exposed 暴露者/潜伏者
    i = 10  # Infected 感染者
    d = 0  # Dead 死亡者
    r = 0  # Recovered 康复者

    alpha = 0.05  # Probability of latent people becoming infected people 潜伏期的人变成感染者的概率
    beta = 0.05  # For susceptible persons, probability of becoming latent after contact with infected persons 对于易感者，与感染者接触后成为潜伏者的可能性
    gamma = 0.08  # Cure rate of infected persons 感染者的治愈率
    delta = 0.02  # Death rate of infected persons 感染者的死亡率
    lamda = 5  # Number of susceptible persons exposed by infected persons per unit time 单位时间内受感染者接触的易受感染人数

    tEnd = 300  # 300天
    t = np.arange(0, tEnd, 1)

    Y0 = (s, e, i, d, r)
    ySEIR = odeint(dySEIDR_C, Y0, t, args=(alpha, beta, gamma, delta, lamda))
    draw_SEIDR(ySEIR)
