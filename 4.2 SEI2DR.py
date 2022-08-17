'''
@coding: utf-8
@author: moyuweiqing

'''

import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def dySEIDR(y, t, alpha1, alpha2, beta, gamma1, gamma2, delta1, delta2, lamda):
    s, e, i1, i2, d, r = y
    n = s + e + i1 + i2 + d + r
    i = i1 + i2

    ds_dt = - lamda * beta * s * (i + e) / n
    de_dt = lamda * beta * s * (i + e) / n - alpha1 * e - alpha2 * e
    di1_dt = alpha1 * e - gamma1 * i1 - delta1 * i1
    di2_dt = alpha2 * e - gamma2 * i2 - delta2 * i2
    dd_dt = delta1 * i1 + delta2 * i2
    dr_dt = gamma1 * i1 + gamma2 * i2

    return np.array([ds_dt, de_dt, di1_dt, di2_dt, dd_dt, dr_dt])

def draw_SEIDR(ySEIR):
    plt.figure(figsize=(16, 9))
    plt.plot(t, ySEIR[:, 0], '--', label='s(t)--SEIR')
    plt.plot(t, ySEIR[:, 1], '-.', label='e(t)--SEIR')
    plt.plot(t, ySEIR[:, 2], '-', label='i1(t)--SEIR')
    plt.plot(t, ySEIR[:, 3], ':', label='i2(t)--SEIR')
    plt.plot(t, ySEIR[:, 4], '.', label='d(t)--SEIR')
    plt.plot(t, ySEIR[:, 5], '--', label='r(t)--SEIR')
    plt.legend()
    # plt.show()
    plt.savefig('./img/' + '4.2 SEI2DR_model.jpg')

if __name__ == '__main__':
    number = 10000      # 模型总人数

    s = 9990            # Susceptible 易感人群数量
    e = 0               # Exposed 暴露者/潜伏者
    i1 = 9              # Infected 感染者
    i2 = 1              # Infected 感染者
    d = 0               # Dead 死亡者
    r = 0               # Recovered 康复者

    alpha1 = 0.04        # Probability of latent people becoming infected people 潜伏期的人变成一般感染者的概率
    alpha2 = 0.01        # Probability of latent people becoming infected people 潜伏期的人变成重症感染者的概率
    beta  = 0.05        # For susceptible persons, probability of becoming latent after contact with infected persons 对于易感者，与感染者接触后成为潜伏者的可能性
    gamma1 = 0.1        # Cure rate of infected persons 一般感染者的治愈率
    gamma2 = 0.02        # Cure rate of infected persons 重症感染者的治愈率
    delta1 = 0.01        # Death rate of infected persons 一般感染者的死亡率
    delta2 = 0.1        # Death rate of infected persons 重症感染者的死亡率
    lamda = 5           # Number of susceptible persons exposed by infected persons per unit time 单位时间内受感染者接触的易受感染人数

    tEnd = 300          # 300天
    t = np.arange(0, tEnd, 1)

    Y0 = (s, e, i1, i2, d, r)
    ySEIR = odeint(dySEIDR, Y0, t, args=(alpha1, alpha2, beta, gamma1, gamma2, delta1, delta2, lamda))
    draw_SEIDR(ySEIR)