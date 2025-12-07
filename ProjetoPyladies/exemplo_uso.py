from dna_sequence import DNASequence
# from fasta_handler import FastaHandler # Comentado para rodar sem arquivo FASTA

def exemplo_basico():
    print("=== EXEMPLO BÁSICO ===")
    dna = DNASequence("ATGCATGCTAGC", "seq1", "Minha sequência de teste")
    print(f"Sequência: {dna}")
    print(f"Comprimento: {len(dna)} bases")
    print(f"GC%: {dna.gc_content():.2f}%\n")


def exemplo_complementaridade():
    print("=== COMPLEMENTARIDADE ===")
    dna = DNASequence("ATGC")
    print(f"Original: {dna.sequence}")
    print(f"Complemento: {dna.complement().sequence}")
    print(f"Reverso complementar: {dna.reverse_complement().sequence}\n")


def exemplo_busca_padroes():
    print("=== BUSCA DE PADRÕES ===")
    dna = DNASequence("ATGATGCATGATG")
    print(f"Sequência: {dna.sequence}")
    print(f"Posições de ATG: {dna.find_all_occurrences('ATG')}")
    print(f"Número de ATGs: {dna.count_pattern_occurrences('ATG')}\n")


def exemplo_analise_completa():
    print("=== ANÁLISE COMPLETA ===")
    # Sequência de comprimento > 14 para usar a fórmula avançada de Tm
    dna = DNASequence("ATGCATGCTAGCTAGCATGCG")
    print(f"ID: {dna.id}")
    print(f"Comprimento: {len(dna)} bases")
    print(f"Composição: {dna.get_base_composition()}")
    print(f"GC%: {dna.gc_content():.2f}%")
    print(f"Tm (Temperatura de Fusão): {dna.calculate_melting_temp():.2f}°C\n")


def exemplo_restricao():
    print("=== SÍTIOS DE RESTRIÇÃO ===")
    dna = DNASequence("ATGAATTCGGCCATGAATTC")
    sites = dna.find_restriction_sites("EcoRI")
    print(f"Sítios EcoRI (GAATTC): {sites}\n")


if __name__ == "__main__":
    exemplo_basico()
    exemplo_complementaridade()
    exemplo_busca_padroes()
    exemplo_analise_completa()
    exemplo_restricao()
    
    print("=== INSTRUÇÕES PARA TESTES ===")
    print("Para rodar os testes de unidade, execute no terminal do VS Code:")
    print("python -m unittest test_dna_sequence.py")