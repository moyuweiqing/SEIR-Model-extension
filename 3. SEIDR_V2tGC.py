'''
@coding: utf-8
@author: moyuweiqing
加入疫苗影响
'''

import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def dySEIDR_V2(y, t, alpha, beta0, beta1, beta2, gamma0, gamma1, gamma2, delta0, delta1, delta2, lamda, v1, v2):
    s0, s1, s2, e0, e1, e2, i0, i1, i2, d0, d1, d2, r = y
    n = s0 + s1 + s2 + e0 + e1 + e2 + i0 + i1 + i2 + d0 + d1 + d2 + r

    if t < 20:             # 20天后开始接种第一针
        v1 = 0
        v2 = 0
    elif t < 50:           # 50天后开始接种第二针
        v2 = 0

    if t >= 30:         # 从第30天开始限制人员接触
        lamda = 2

    if t >= 30 and t <= 60:                   # 从第30天开始病毒发生变异
        beta0 = beta0 * 1.2                   # 传染性增强
        beta1 = beta1 * 1.2                   # 传染性增强
        beta2 = beta2 * 1.2                   # 传染性增强
        delta0 = delta0 * 0.8                 # 致死性降低
        delta1 = delta1 * 0.8                 # 致死性降低
        delta2 = delta2 * 0.8                 # 致死性降低
    elif t > 60:                              # 第二次变异
        beta0 = beta0 * 1.5                   # 传染性进一步增强
        beta1 = beta1 * 1.5                   # 传染性进一步增强
        beta2 = beta2 * 1.5                   # 传染性进一步增强
        delta0 = delta0 * 0.5                 # 致死性进一步降低
        delta1 = delta1 * 0.5                 # 致死性进一步降低
        delta2 = delta2 * 0.5                 # 致死性进一步降低

    i_e = i0 + i1 + i2 + e0 + e1 + e2

    ds0_dt = - lamda * beta0 * s0 * i_e / n - v1 * s0
    ds1_dt = - lamda * beta1 * s1 * i_e / n + v1 * s0 - v2 * s1
    ds2_dt = - lamda * beta2 * s2 * i_e / n + v2 * s1
    de0_dt = lamda * beta0 * s0 * i_e / n - alpha * e0
    de1_dt = lamda * beta1 * s1 * i_e / n - alpha * e1
    de2_dt = lamda * beta2 * s2 * i_e / n - alpha * e2
    di0_dt = alpha * e0 - gamma0 * i0 - delta0 * i0
    di1_dt = alpha * e1 - gamma1 * i1 - delta1 * i1
    di2_dt = alpha * e2 - gamma2 * i2 - delta2 * i2
    dd0_dt = delta0 * i0
    dd1_dt = delta1 * i1
    dd2_dt = delta2 * i2
    dr_dt = gamma0 * i0 + gamma1 * i1 + gamma2 * i2

    return np.array([ds0_dt, ds1_dt, ds2_dt, de0_dt, de1_dt, de2_dt, di0_dt, di1_dt, di2_dt, dd0_dt, dd1_dt, dd2_dt, dr_dt])

def draw_SEIDR(ySEIR):
    plt.figure(figsize=(16, 9))
    plt.plot(t, ySEIR[:, 0] + ySEIR[:, 1] + ySEIR[:, 2], '--', label='s(t)--SEIR')
    plt.plot(t, ySEIR[:, 3] + ySEIR[:, 4] + ySEIR[:, 5], '-.', label='e(t)--SEIR')
    plt.plot(t, ySEIR[:, 6] + ySEIR[:, 7] + ySEIR[:, 8], '-', label='i(t)--SEIR')
    plt.plot(t, ySEIR[:, 9] + ySEIR[:, 10] + ySEIR[:, 11], '.', label='d(t)--SEIR')
    plt.plot(t, ySEIR[:, 12], '--', label='r(t)--SEIR')
    plt.legend()
    # plt.show()
    plt.savefig('./img/' + '3. SEIDR_V2tGC_model.jpg')

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
    lamda = 5           # Number of susceptible persons exposed by infected persons per unit time 单位时间内受感染者接触的易受感染人数

    tEnd = 300          # 300天
    t = np.arange(0, tEnd, 1)


    v1 = 0.02            # Vaccination coverage rates of first dose 第一针疫苗接种率
    v2 = 0.02            # Vaccination coverage rates of second dose 第二针疫苗接种率
    beta0 = 0.05        # The probability of becoming latent in unvaccinated population 未打疫苗人群的患病概率
    beta1 = 0.01        # The probability of becoming latent in vaccinated population 已打疫苗人群的感染率
    beta2 = 0.005       # The probability of becoming latent in vaccinated population 已打疫苗人群的感染率
    gamma0 = 0.08       # Cure rates of unvaccinated infected persons 未接种疫苗的感染者的治愈率
    gamma1 = 0.10       # The cure rate of vaccinated infected people 接种了疫苗的感染者的治愈率
    gamma2 = 0.10       # The cure rate of vaccinated infected people 接种了疫苗的感染者的治愈率
    delta0 = 0.02       # Death rate of unvaccinated infected persons 未接种疫苗的感染者的死亡率
    delta1 = 0.001      # Death rate of vaccinated infected persons 接种了疫苗的感染者的死亡率
    delta2 = 0.0008     # Death rate of vaccinated infected persons 接种了疫苗的感染者的死亡率
    s0 = 9990           # Number of unvaccinated people 未接种疫苗人群
    s1 = 0              # The number of people who receive one dose of vaccine 接种一针疫苗人群
    s2 = 0              # The number of people who receive two doses of vaccine 接种两针疫苗人群
    e0 = 0              # The latent people who do not receive vaccine 未接种疫苗的潜伏人群
    e1 = 0              # The latent people who receive one dose of vaccine 接种一针疫苗的潜伏人群
    e2 = 0              # The latent people who receive two doses of vaccine 接种两针疫苗的潜伏人群
    i0 = 10             # The infected people who do not receive vaccine 未接种疫苗的感染人群
    i1 = 0              # The infected people who receive one dose of vaccine 接种一针疫苗的感染人群
    i2 = 0              # The infected people who receive two doses of vaccine 接种两针疫苗的感染人群
    d0 = 0              # The dead people who do not receive vaccine 未接种疫苗的死亡人群
    d1 = 0              # The dead people who receive one dose of vaccine 接种一针疫苗的死亡人群
    d2 = 0              # The dead people who receive two doses of vaccine 接种两针疫苗的感染人群

    Y0 = (s0, s1, s2, e0, e1, e2, i0, i1, i2, d0, d1, d2, r)
    ySEIR = odeint(dySEIDR_V2, Y0, t, args=(alpha, beta0, beta1, beta2, gamma0, gamma1, gamma2, delta0, delta1, delta2, lamda, v1, v2))
    draw_SEIDR(ySEIR)