import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.fft import fft, fftfreq

import os


def load_data(file_path):
    """
    Carrega os dados do arquivo e retorna um DataFrame pandas

    Args:
        file_path (str): Caminho para o arquivo de dados

    Returns:
        pandas.DataFrame: DataFrame contendo os dados carregados
    """
    # Carregar os dados do arquivo
    try:
        # Definir nomes das colunas
        column_names = ['Reynolds', 'Cd', 'Cl']

        # Carregar dados sem cabeçalho
        df = pd.read_csv(file_path, sep='\t', header=None, names=column_names)

        print(f"Dados carregados com sucesso. Formato: {df.shape}")
        return df
    except Exception as e:
        print(f"Erro ao carregar os dados: {e}")
        return None

def plot_time_series_by_reynolds(df, output_dir='plots',cut=0):
    """
    Cria gráficos de série temporal para Cd e Cl separados por Reynolds,
    incluindo linhas de média para ambos os coeficientes

    Args:
        df (pandas.DataFrame): DataFrame contendo os dados
        output_dir (str): Diretório para salvar os gráficos
    """
    # Criar diretório para os gráficos se não existir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    #Vetor temporario
    dados_cd = []

    # Obter valores únicos de Reynolds
    reynolds_values = df['Reynolds'].unique()

    # Criar gráficos de série temporal para cada Reynolds
    for reynolds in reynolds_values:
        # Filtrar dados para o Reynolds atual
        df_reynolds = df[df['Reynolds'] == reynolds]

        if len(df_reynolds) <= cut:
            print(f"-->Aviso: Reynolds {reynolds} tem menos de {cut} dados. Pulando.")
            continue

        # Adicionar coluna de índice para representar a sequência temporal
        df_temp = df_reynolds.copy()
        df_temp['Tempo'] = df_temp.index

        # Calcular médias (após o ponto 2000)
        mean_cd = np.mean(df_temp['Cd'][cut:])
        mean_cl = np.mean(df_temp['Cl'][cut:])

        # Gráfico de série temporal para Cd
        plt.figure(figsize=(12, 6))
        plt.plot(df_temp['Tempo'][cut:], df_temp['Cd'][cut:], 'b-', linewidth=1, label='Cd')
        plt.axhline(y=mean_cd, color='red', linestyle='--', label=f'Média (Cd = {mean_cd:.4f})')
        plt.title(f'Variação de Cd ao longo do tempo para Reynolds = {reynolds}')
        plt.xlabel('Sequência')
        plt.ylabel('Cd')
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'{output_dir}/time_series_cd_reynolds_{int(reynolds)}.png', dpi=300)
        plt.close()

        # Gráfico de série temporal para Cl
        plt.figure(figsize=(12, 6))
        plt.plot(df_temp['Tempo'][cut:], df_temp['Cl'][cut:], 'g-', linewidth=1, label='Cl')
        plt.axhline(y=mean_cl, color='orange', linestyle='--', label=f'Média (Cl = {mean_cl:.4f})')
        plt.title(f'Variação de Cl ao longo do tempo para Reynolds = {reynolds}')
        plt.xlabel('Sequência')
        plt.ylabel('Cl')
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'{output_dir}/time_series_cl_reynolds_{int(reynolds)}.png', dpi=300)
        plt.close()

        dados_cd.append((int(reynolds),mean_cd))

    print(f"Gráficos de série temporal por Reynolds salvos no diretório '{output_dir}'")
    return dados_cd

def analyze_fft_cl_by_reynolds(df, timestep=1.0, cut=0, D=40, U=0.0557, output_dir='plots'):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    reynolds_values = df['Reynolds'].unique()
    dados_st = []

    for reynolds in reynolds_values:
        df_re = df[df['Reynolds'] == reynolds]
        cl_series = df_re['Cl'].values

        if len(cl_series) <= cut:
            print(f"-->Aviso: Reynolds {reynolds} tem menos de {cut} dados. Pulando.")
            continue

        cl_series_cut = cl_series[cut:]
        n = len(cl_series_cut)

        window = np.hanning(n)
        cl_windowed = cl_series_cut * window

        yf = fft(cl_windowed)
        xf = fftfreq(n, d=timestep)[:n//2]
        spectrum = 2.0/n * np.abs(yf[:n//2])

        peak_index = np.argmax(spectrum[1:]) + 1
        freq_dominante = xf[peak_index]
        strouhal = (freq_dominante * D) / U

        plt.figure(figsize=(10, 6))
        plt.plot(xf, spectrum)
        plt.title(f'FFT de Cl - Re = {reynolds} | f = {freq_dominante:.4f} | St = {strouhal:.4f}')
        plt.xlabel('Frequência (1/tempo)')
        plt.ylabel('Amplitude')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(f"{output_dir}/fft_cl_reynolds_{int(reynolds)}.png", dpi=300)
        plt.close()

        print(f"Re = {reynolds}: f = {freq_dominante:.6f}, St = {strouhal:.6f}")
        dados_st.append((int(reynolds), strouhal))

    return dados_st

def main():
    file_path = f"/home/thomas/remote/hal/Re300/EmerichCube/force.txt"
    dir_name = "300Emerich"
    cut = 0

    if not os.path.exists(file_path):
        print(f"Erro: O arquivo '{file_path}' não foi encontrado.")
        return

    df = load_data(file_path)
    if df is None:
        return

    dados_cd = plot_time_series_by_reynolds(df, output_dir= dir_name,cut=cut)
    dados_st = analyze_fft_cl_by_reynolds(df, output_dir= dir_name,timestep=100, cut=cut, D=40, U=0.0577)

    # Combinar dados por Reynolds
    cd_dict = dict(dados_cd)
    st_dict = dict(dados_st)

    dados_finais = []
    for Re in sorted(cd_dict.keys()):
        if Re in st_dict:
            Cd = cd_dict[Re]
            St = st_dict[Re]
            dados_finais.append((Re, Cd, St))

    structured_array = np.array(dados_finais, dtype=[("Re", "i4"), ("Cd", "f4"), ("St", "f4")])
    np.save(f"{dir_name}/Cd_St_Re.npy", structured_array)

    print("Dados salvos com sucesso em 'Cd_St_Re.npy'")

if __name__ == "__main__":
    main()
