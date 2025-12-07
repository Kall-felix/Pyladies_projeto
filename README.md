# 🧬 DNA Toolkit: Análise de Sequências em Python

Este projeto é um *toolkit* de bioinformática implementado em Python puro, focado na modelagem, manipulação e análise de sequências de DNA. Ele foi desenvolvido para aplicar boas práticas de programação (como o princípio da Responsabilidade Única - SRP) e testabilidade, conforme as diretrizes do nosso desafio.

### 🔬 Funcionalidades Chave

O pacote oferece as seguintes capacidades analíticas básicas e avançadas:

* **Composição de Bases:** Contagem de bases (A, T, G, C, N).
* **Conteúdo GC:** Cálculo do percentual de Guanina e Citosina.
* **Complementaridade:** Geração do Complemento Reverso da sequência, crucial em biologia molecular.
* **Transcrição:** Conversão de DNA para RNA (T $\rightarrow$ U).
* **Termostabilidade (Tm):** Estimativa da Temperatura de Fusão (Tm), essencial para *design* de *primers*.
* **Sítios de Restrição:** Localização de sítios de enzimas comuns (e.g., EcoRI).
* **Busca de Padrões:** Encontra todas as ocorrências de um padrão específico na sequência.
* **Manipulação de Arquivos:** Leitura e escrita de sequências no formato padrão FASTA.

---

## 🛠️ Estrutura e Modelagem do Projeto

O projeto segue um design modular, aplicando princípios de boa arquitetura para garantir que cada componente tenha uma responsabilidade clara (**SRP**):

* **`dna_sequence.py`:** Contém a classe principal **`DNASequence`** e o utilitário **`DNAValidator`**. A responsabilidade é exclusiva sobre a lógica, manipulação e análise da sequência.
* **`fasta_handler.py`:** Contém a classe **`FastaHandler`**. A responsabilidade é estritamente sobre a entrada e saída de dados de arquivos FASTA.
* **`test_dna_sequence.py`:** Contém a suíte completa de testes de unidade.
* **`exemplo_uso.py`:** Script de demonstração das funcionalidades.

---

## 🚀 Como Rodar e Testar

### Pré-requisitos

O projeto usa apenas bibliotecas nativas do Python. Você só precisa ter o **Python 3.x** instalado.

### 1. Demonstração de Uso

Para ver todos os métodos analíticos em ação, execute o script de exemplo:

```bash
python exemplo_uso.py
