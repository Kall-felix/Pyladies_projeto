from typing import List
from pathlib import Path


class FileFormatError(Exception):
    pass


class FastaHandler:
    """Manipulação de arquivos FASTA"""
    
    @staticmethod
    def read(filepath: str, seq_class):
        """Lê arquivo FASTA e retorna lista de sequências."""
        if not Path(filepath).exists():
            raise FileFormatError(f"Arquivo não encontrado: {filepath}")
        
        sequences = []
        current_id = None
        current_desc = None
        current_seq_parts = []
        
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                if line.startswith('>'):
                    # Se houver uma sequência anterior, processe-a antes de iniciar a nova
                    if current_id:
                        seq = ''.join(current_seq_parts)
                        sequences.append(seq_class(seq, current_id, current_desc))
                    
                    # Inicia nova sequência
                    header = line[1:].split(maxsplit=1)
                    current_id = header[0]
                    current_desc = header[1] if len(header) > 1 else ""
                    current_seq_parts = []
                else:
                    current_seq_parts.append(line)
            
            # Processa a última sequência no arquivo
            if current_id:
                seq = ''.join(current_seq_parts)
                sequences.append(seq_class(seq, current_id, current_desc))
        
        return sequences
    
    @staticmethod
    def write(sequences: List, filepath: str, width: int = 80):
        """Escreve sequências em arquivo FASTA."""
        with open(filepath, 'w') as f:
            for seq in sequences:
                f.write(f">{seq.id}")
                if seq.description:
                    f.write(f" {seq.description}")
                f.write("\n")
                
                sequence = seq.sequence
                for i in range(0, len(sequence), width):
                    f.write(sequence[i:i+width] + "\n")


if __name__ == "__main__":
    from dna_sequence import DNASequence
    
    # Crie um arquivo 'sequences.fasta' de exemplo para este teste
    print("Para testar este bloco, crie um arquivo 'sequences.fasta' no diretório.")
    
    # Exemplo: Ler arquivo
    # try:
    #     sequences = FastaHandler.read("sequences.fasta", DNASequence)
    #     for seq in sequences:
    #         print(f"{seq.id}: GC% = {seq.gc_content():.2f}%")
    # except FileFormatError as e:
    #     print(e)