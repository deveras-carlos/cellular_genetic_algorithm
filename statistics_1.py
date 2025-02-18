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
    funcoes = sorted(set(res['i'] for res in resultados))

    for func in funcoes:
        func_resultados = [res for res in resultados if res['i'] == func]
        
        # Generate plots based on the number of variables
        variaveis_set = sorted(set(res['variaveis'] for res in func_resultados))
        for variaveis in variaveis_set:
            var_resultados = [res for res in func_resultados if res['variaveis'] == variaveis]
            
            plt.figure(figsize=(10, 5))
            tentativas_set = sorted(set(res['tentativas'] for res in var_resultados))
            for tentativas in tentativas_set:
                trial_resultados = [res for res in var_resultados if res['tentativas'] == tentativas]
                if trial_resultados:
                    pops = [res['populacao'] for res in trial_resultados]
                    tempos = [res['tempo_execucao'] for res in trial_resultados]
                    plt.plot(pops, tempos, marker='o', label=f'Tentativas: {tentativas}')
            
            plt.xlabel("Tamanho da População")
            plt.ylabel("Tempo de Execução (s)")
            plt.title(f"Função {func} - Variáveis: {variaveis}")
            plt.legend()
            output_filename = f"tempos_execucao_func_{func}_vars_{variaveis}.png"
            plt.savefig(output_filename)
            plt.close()

        # Generate plots based on the number of trials
        tentativas_set = sorted(set(res['tentativas'] for res in func_resultados))
        for tentativas in tentativas_set:
            trial_resultados = [res for res in func_resultados if res['tentativas'] == tentativas]
            
            plt.figure(figsize=(10, 5))
            variaveis_set = sorted(set(res['variaveis'] for res in trial_resultados))
            for variaveis in variaveis_set:
                var_resultados = [res for res in trial_resultados if res['variaveis'] == variaveis]
                if var_resultados:
                    pops = [res['populacao'] for res in var_resultados]
                    tempos = [res['tempo_execucao'] for res in var_resultados]
                    plt.plot(pops, tempos, marker='o', label=f'Variáveis: {variaveis}')
            
            plt.xlabel("Tamanho da População")
            plt.ylabel("Tempo de Execução (s)")
            plt.title(f"Função {func} - Tentativas: {tentativas}")
            plt.legend()
            output_filename = f"tempos_execucao_func_{func}_trials_{tentativas}.png"
            plt.savefig(output_filename)
            plt.close()

def main():
    resultados = processar_resultados()
    salvar_csv(resultados)
    plotar_tempos_execucao(resultados)
    print(f"Resultados salvos em {OUTPUT_CSV} e gráfico de tempo de execução gerado.")

if __name__ == "__main__":
    main()
