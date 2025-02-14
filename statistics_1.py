import os
import csv
import re
import matplotlib.pyplot as plt

RESULTS_DIR = "results"
OUTPUT_CSV = "resumo.csv"

def parse_file(filepath):
    data = {}
    with open(filepath, 'r') as file:
        content = file.read()
        
        # Extrair dados necessários
        data['tempo_execucao'] = round(float(re.search(r'Execution Time: ([\d\.]+)', content).group(1)), 4)
        data['tamanho_populacao'] = int(re.search(r'Population Size: (\d+)', content).group(1))
        data['tamanho_individuo'] = int(re.search(r'Individual Size: (\d+)', content).group(1))
        data['melhor_fitness'] = round(float(re.search(r'Best Fitness: ([\d\.]+)', content).group(1)), 4)
        data['melhor_geracao'] = int(re.search(r'Best Generation: (\d+)', content).group(1))
        
        # Extrair valores do nome do arquivo
        match = re.search(r'results_i(\d+)_vars_(\d+)_trials(\d+)_pop(\d+)_j(\d+)\.txt', filepath)
        if match:
            data['i'] = int(match.group(1))
            data['variaveis'] = int(match.group(2))
            data['tentativas'] = int(match.group(3))
            data['populacao'] = int(match.group(4))
            data['j'] = int(match.group(5))
    
    return data

def processar_resultados():
    melhores_resultados = {}
    
    for filename in os.listdir(RESULTS_DIR):
        if filename.endswith(".txt"):
            filepath = os.path.join(RESULTS_DIR, filename)
            dados = parse_file(filepath)
            
            chave = (dados['i'], dados['variaveis'], dados['tentativas'], dados['populacao'])
            
            if chave not in melhores_resultados or dados['melhor_fitness'] < melhores_resultados[chave]['melhor_fitness']:
                melhores_resultados[chave] = dados
    
    return list(melhores_resultados.values())

def salvar_csv(resultados):
    with open(OUTPUT_CSV, 'w', newline='') as csvfile:
        cabecalhos = ["i", "variaveis", "tentativas", "populacao", "melhor_fitness", "melhor_geracao", "tempo_execucao"]
        escritor = csv.DictWriter(csvfile, fieldnames=cabecalhos)
        escritor.writeheader()
        for linha in resultados:
            linha_filtrada = {chave: linha[chave] for chave in cabecalhos}
            escritor.writerow(linha_filtrada)

def plotar_tempos_execucao(resultados):
    tentativas = [res['tentativas'] for res in resultados]
    tempos = [res['tempo_execucao'] for res in resultados]
    
    plt.figure(figsize=(10, 5))
    plt.bar(tentativas, tempos, color='blue')
    plt.xlabel("Número de Tentativas")
    plt.ylabel("Tempo de Execução (s)")
    plt.title("Tempo de Execução por Número de Tentativas")
    plt.xticks(tentativas)
    plt.savefig("tempos_execucao.png")
    plt.show()

def main():
    resultados = processar_resultados()
    salvar_csv(resultados)
    plotar_tempos_execucao(resultados)
    print(f"Resultados salvos em {OUTPUT_CSV} e gráfico de tempo de execução gerado.")

if __name__ == "__main__":
    main()
