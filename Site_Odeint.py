import streamlit as st
from scipy.integrate import odeint
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.animation import FuncAnimation
from seaborn.palettes import blend_palette

# Nome e descrição
st.title('Simulação de um trocador de calor tubular simples')
st.write('Este site simula um trocador de calor tubular simples que aquece um fluido conforme o mesmo passa por ele. Ao rodar a simulação, você poderá visualizar o perfil de temperatura do fluido no trocador de calor conforme o tempo passa. Você também poderá visualizar o gráfico de variação da temperatura quando o trocador atinge o regime permanente.')

# Carregar a imagem
st.image('Caso 2.png', use_column_width=True)
st.write('Figura exemplificando o trocador. Autoria própria.')

# Valores input
L = st.number_input('Comprimento do tubo (m)', min_value=0)
r = st.number_input('Raio do tubo (m)', min_value=0.0)
n = st.number_input('Número de nós para a discritização', min_value=1)
m = st.number_input('Vazão Mássica (kg/s)', min_value=0.0)
Cp = st.number_input('Capacidade de calor específico do fluido (J/kg.K)', min_value=0.0)
rho = st.number_input('Massa específica do fluido (kg/m³)', min_value=0.0)
Ti = st.number_input('Temperatura inicial do tubo (K)')
T0 = st.number_input('Temperatura inicial do fluido (K)')
q_fluxo = st.number_input('Fluxo de calor (W/m²)', min_value=0.0)
t_final = st.number_input('Tempo de simulação (s)', min_value=0.0)
dt = st.number_input('Passo de tempo (s)', min_value=0.0)

# Calculo de dx
dx = L / n

# Criação de um botão para rodar a simulação
if st.button('Rodar Simulação'):
    x = np.linspace(dx/2, L-dx/2, n)
    T = np.ones(n) * T0
    dTdt = np.zeros(n)
    t = np.arange(0, t_final, dt)
    T_out = [T]

    # Cálculo de T em regime permanente
    T_permanente = Ti + q_fluxo * 2 * np.pi * r * x / (m * Cp)
    
    # Criando a figura para o gráfico em regime permanente
    fig_permanente = plt.figure(figsize=(5, 4))
    
    #Criando figura para a animação
    fig_animacao = plt.figure(figsize=(5, 4))
    
    # Função que define a EDO para a variação da temperatura
    def dTdt_function(T, t):
        dTdt = np.zeros(n)
        dTdt[1:n] = (m*Cp*(T[0:n-1]-T[1:n])+q_fluxo*2*np.pi*r*dx)/(rho*Cp*dx*np.pi*r**2)
        dTdt[0] = (m*Cp*(Ti-T[0])+q_fluxo*2*np.pi*r*dx)/(rho*Cp*dx*np.pi*r**2)
        return dTdt

    # Resolvendo a EDO usando o odeint
    T_out = odeint(dTdt_function, T, t)
    T_out = T_out - 273.15
    
    # Criação do DataFrame
    df_Temp = pd.DataFrame(np.array(T_out), columns=x)

    # Create uma paleta de cores
    paleta_calor = blend_palette(['yellow', 'orange','red'], as_cmap=True, n_colors=100)

    # Função que atualiza o plot
    def update_plot(t):
        plt.clf()
        line = pd.DataFrame(df_Temp.iloc[t, :]).T
        sns.heatmap(line, cmap=paleta_calor)
        plt.title(f'Tempo: {t} (s)')

    # Criando a animação
    ani = FuncAnimation(fig_animacao, update_plot, frames=df_Temp.shape[0], repeat=False)

    # Salvar a animação como gif
    ani.save('Variação da Temperatura - Caso 2.gif', writer='pillow', fps=10)
    
    # Exibindo a simulação
    with st.expander("Visualização da Simulação em tempo real (Clique aqui para ver)"):
        st.image('Variação da Temperatura - Caso 2.gif')
    
    #Exibindo o gráfico de variação da temperatura ao longo do comprimento em regime permanente
    st.write('Gráfico da variação da temperatura ao longo do comprimento em regime permanente')
    plt.figure(fig_permanente)
    plt.plot(x, T_permanente, color='blue')  
    plt.xlabel('Comprimento (x)')
    plt.ylabel('Temperatura (°C)')
    st.pyplot(plt)






