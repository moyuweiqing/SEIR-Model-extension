'''
@coding: utf-8
@author: moyuweiqing

加入疫苗影响
'''

import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def dySEIDR_Vt(y, t, alpha, beta0, beta1, gamma0, gamma1, delta0, delta1, lamda, v):
    s0, s1, e0, e1, i0, i1, d0, d1, r = y
    n = s0 + s1 + e0 + e1 + i0 + i1 + d0 + d1 + r

    if t < 20:             # 20天后开始接种
        v = 0

    ds0_dt = - lamda * beta0 * s0 * ( i0 + i1 + e0 + e1 ) / n - v * s0
    ds1_dt = - lamda * beta1 * s1 * ( i0 + i1 + e0 + e1 ) / n + v * s0
    de0_dt = lamda * beta0 * s0 * ( i0 + i1 + e0 + e1 ) / n - alpha * e0
    de1_dt = lamda * beta1 * s1 * ( i0 + i1 + e0 + e1 ) / n - alpha * e1
    di0_dt = alpha * e0 - gamma0 * i0 - delta0 * i0
    di1_dt = alpha * e1 - gamma1 * i1 - delta1 * i1
    dd0_dt = delta0 * i0
    dd1_dt = delta1 * i1
    dr_dt = gamma0 * i0 + gamma1 * i1

    return np.array([ds0_dt, ds1_dt, de0_dt, de1_dt, di0_dt, di1_dt, dd0_dt, dd1_dt, dr_dt])

def draw_SEIDR(ySEIR):
    plt.figure(figsize=(16, 9))
    plt.plot(t, ySEIR[:, 0] + ySEIR[:, 1], '--', label='s(t)--SEIR')
    plt.plot(t, ySEIR[:, 2] + ySEIR[:, 3], '-.', label='e(t)--SEIR')
    plt.plot(t, ySEIR[:, 4] + ySEIR[:, 5], '-', label='i(t)--SEIR')
    plt.plot(t, ySEIR[:, 6] + ySEIR[:, 7], '.', label='d(t)--SEIR')
    plt.plot(t, ySEIR[:, 8], '--', label='r(t)--SEIR')
    plt.legend()
    # plt.show()
    plt.savefig('./img/' + '2.1.2 SEIDR_Vt_model.jpg')

if __name__ == '__main__':
    number = 10000      # 模型总人数

    s = 9990            # Susceptible 易感人群数量
    e = 0               # Exposed 暴露者/潜伏者
    i = 10              # Infected 感染者
    d = 0               # Dead 死亡者
    r = 0               # Recovered 康复者

    alpha = 0.05        # Probability of latent people becoming infected people 潜伏期的人变成感染者的概率
    beta  = 0.05        # For susceptible persons, probability of becoming latent after contact with infected persons 对于易感者，与感染者接触后成为潜伏者的可能性
    gamma = 0.08        # Cure rate of infected persons 感染者的治愈率
    delta = 0.02        # Death rate of infected persons 感染者的死亡率
    lamda = 5          # Number of susceptible persons exposed by infected persons per unit time 单位时间内受感染者接触的易受感染人数

    tEnd = 300          # 300天
    t = np.arange(0, tEnd, 1)


    v = 0.02            # Vaccination coverage rates 疫苗接种率
    beta0 = 0.05        # The probability of becoming latent in unvaccinated population 未打疫苗人群的患病概率
    beta1 = 0.01        # The probability of becoming latent in vaccinated population 已打疫苗人群的感染率
    gamma0 = 0.08       # Cure rates of unvaccinated infected persons 未接种疫苗的感染者的治愈率
    gamma1 = 0.10       # The cure rate of vaccinated infected people 接种了疫苗的感染者的治愈率
    delta0 = 0.02       # Death rate of unvaccinated infected persons 未接种疫苗的感染者的死亡率
    delta1 = 0.001      # Death rate of vaccinated infected persons 接种了疫苗的感染者的死亡率
    s0 = 9990
    s1 = 0
    e0 = 0
    e1 = 0
    i0 = 10
    i1 = 0
    d0 = 0
    d1 = 0

    Y0 = (s0, s1, e0, e1, i0, i1, d0, d1, r)
    ySEIR = odeint(dySEIDR_Vt, Y0, t, args=(alpha, beta0, beta1, gamma0, gamma1, delta0, delta1, lamda, v))
    draw_SEIDR(ySEIR)