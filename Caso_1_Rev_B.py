import streamlit as st
from scipy.integrate import odeint
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.animation import FuncAnimation
from seaborn.palettes import blend_palette


def run_simulation(L, r, n, m, Cp, rho, Ti, T0, q_fluxo, t_final, dt):
    dx = L / n
    x = np.linspace(dx/2, L-dx/2, n)
    T = np.ones(n) * T0
    t = np.arange(0, t_final, dt)
    
    # Criando a figura para o gráfico em regime permanente
    fig_permanente = plt.figure(figsize=(8, 6))

    # Função que define a EDO para a variação da temperatura
    def dTdt_function(T, t):
        dTdt = np.zeros(n)
        dTdt[1:n] = (m*Cp*(T[0:n-1]-T[1:n])+q_fluxo*2*np.pi*r*dx)/(rho*Cp*dx*np.pi*r**2)
        dTdt[0] = (m*Cp*(Ti-T[0])+q_fluxo*2*np.pi*r*dx)/(rho*Cp*dx*np.pi*r**2)
        return dTdt
    
    # Resolvendo a EDO usando o odeint
    T_out = odeint(dTdt_function, T, t)
    T_out = T_out
    
    # Criação do DataFrame
    df_Temp = pd.DataFrame(np.array(T_out), columns=x)
    
    # Criando uma paleta de cores
    paleta_calor = blend_palette(['yellow', 'orange','red'], as_cmap=True, n_colors=100)
    
    # Função que atualiza o plot
    def update_plot(t):
        plt.clf()
        line = pd.DataFrame(df_Temp.iloc[t, :]).T
        sns.heatmap(line, cmap=paleta_calor)
        plt.title(f'Tempo: {t} (s)')

    #Criando figura para a animação
    fig_animacao = plt.figure(figsize=(8, 6))
    
    # Criando a animação
    ani = FuncAnimation(fig_animacao, update_plot, frames=df_Temp.shape[0], repeat=False)

    # Salvar a animação como gif
    save = ani.save('Variação da Temperatura - Caso I.gif', writer='pillow', fps=10)
    
    # Exibindo a simulação
    with st.expander("Visualização da Simulação em tempo real para o fluido (Clique aqui para ver)"):
        st.write('Variação da temperatura do fluido que passa pelo trocador ao longo do tempo e do comprimento.')
        st.write('Tempo representado acima do GIF, em segundos. Temperaturas em Kelvin representadas na escala variável do eixo y. Comprimento do trocador representado no eixo x do GIF.')
        st.image('Variação da Temperatura - Caso I.gif')

    #Exibindo o gráfico de variação da temperatura ao longo do comprimento em regime permanente
    plt.figure(fig_permanente)
    plt.plot(x, df_Temp.iloc[-1, :], color='blue')  
    plt.xlabel('Comprimento (m)')
    plt.ylabel('Temperatura (K)')
    st.pyplot(plt)

st.title('Simulador TROCAL - Simulação de um trocador de calor tubular simples')
st.write('Este é um simulador de um trocador de calor tubular simples que aquece um fluido conforme o mesmo passa por ele. Ao rodar a simulação, você poderá visualizar o perfil de temperatura do fluido ao longo do trocador conforme o tempo passa. Você também poderá visualizar as temperaturas em regime permanente ao longo do comprimento do trocador.')
st.write('Este simulador utiliza a seguinte equação de balanço de energia:')
st.image('Equacao Caso 1.png', use_column_width=True)
st.write('ATENÇÃO: Ao final desta página, você também encontrará um botão que roda a simulação com um exemplo pré-definido ("Rodar exemplo padrão"). Este exemplo leva em torno de 3 minutos para rodar. Caso queira usar seus próprios valores, use o botão "Rodar simulação" e recomenda-se utilizar um número de nós pelo menos igual ou maior que 10, para melhor visualização dos resultados.')

# Carregar a imagem
st.image('Caso 1.png', use_column_width=True)
st.write('Figura exemplificando o trocador. Autoria própria.')

# Valores input
L = st.number_input('Comprimento do tubo (m)', min_value=0)
r = st.number_input('Raio do tubo (m)', min_value=0.0)
n = st.number_input('Número de nós para a discritização', min_value=1)
m = st.number_input('Vazão Mássica (kg/s)', min_value=0.0)
Cp = st.number_input('Capacidade de calor específico do fluido (J/kg.K)', min_value=0.0)
rho = st.number_input('Massa específica do fluido (kg/m³)', min_value=0.0)
Ti = st.number_input('Temperatura de entrada do fluido (K)')
T0 = st.number_input('Temperatura inicial do trocador (K)')
q_fluxo = st.number_input('Fluxo de calor (W/m²)', min_value=0.0)
t_final = st.number_input('Tempo de simulação (s)', min_value=0.0)
dt = st.number_input('Passo de tempo (s)', min_value=0.0)

if st.button('Rodar Simulação'):
    run_simulation(L, r, n, m, Cp, rho, Ti, T0, q_fluxo, t_final, dt)
elif st.button('Rodar exemplo padrão'):
    run_simulation(25, 0.1, 100, 3, 4180, 1000, 400, 300, 100000, 350, 1)
