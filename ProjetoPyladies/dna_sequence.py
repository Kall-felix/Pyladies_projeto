from typing import Optional, List, Dict, Tuple

# --- Exceções e Validação ---

class InvalidSequenceError(Exception):
    """Exceção levantada para sequências de DNA que contêm bases inválidas."""
    pass


class DNAValidator:
    """Utilitário estático para validar bases de DNA."""
    _VALID_BASES = {'A', 'T', 'G', 'C', 'N'}
    
    @staticmethod
    def is_valid(sequence: str) -> bool:
        """Verifica se a sequência é não-vazia e contém apenas bases válidas."""
        if not sequence:
            return False
        # Uso de set.issubset para verificação eficiente
        return set(sequence.upper()) <= DNAValidator._VALID_BASES


# --- Classe Principal ---

class DNASequence:
    """
    Representa e manipula uma sequência de DNA.

    Oferece métodos para operações básicas (composição de bases, GC content),
    complementaridade (complemento reverso), transcrição e análise.
    """
    
    _COMPLEMENT_MAP = str.maketrans('ATGCN', 'TACGN')
    
    def __init__(self, sequence: str, seq_id: Optional[str] = None, 
                 description: Optional[str] = None):
        """Inicializa a sequência, validando e padronizando-a."""
        self._sequence = sequence.upper().strip()
        self.id = seq_id if seq_id is not None else "unnamed"
        self.description = description if description is not None else ""
        
        if not DNAValidator.is_valid(self._sequence):
            raise InvalidSequenceError("A sequência contém caracteres não-DNA inválidos.")
    
    @property
    def sequence(self) -> str:
        """Retorna a sequência de DNA padronizada."""
        return self._sequence
    
    def __len__(self) -> int:
        return len(self._sequence)
    
    def __str__(self) -> str:
        return f"DNASequence(id='{self.id}', length={len(self)})"
    
    def __repr__(self) -> str:
        preview = self._sequence[:50] + "..." if len(self) > 50 else self._sequence
        return f"DNASequence(id='{self.id}', seq='{preview}')"
    
    def __getitem__(self, key):
        """Permite fatiar a sequência. Retorna uma nova DNASequence se for um slice."""
        result = self._sequence[key]
        if isinstance(key, slice):
            return self._create_slice_sequence(result)
        return result
    
    def __eq__(self, other) -> bool:
        """Compara duas sequências de DNA pelo seu conteúdo."""
        if not isinstance(other, DNASequence):
            return NotImplemented
        return self._sequence == other._sequence
    
    def _create_slice_sequence(self, sub_sequence: str) -> 'DNASequence':
        """Método auxiliar para criar um objeto DNASequence a partir de um slice."""
        new_id = f"{self.id}_slice_{len(sub_sequence)}"
        return DNASequence(sub_sequence, new_id)

    # --- OPERAÇÕES BÁSICAS ---
    
    def get_base_composition(self) -> Dict[str, int]:
        """Calcula a contagem de cada base (A, T, G, C, N) na sequência."""
        bases = DNAValidator._VALID_BASES
        return {base: self._sequence.count(base) for base in sorted(list(bases))}
    
    def gc_content(self) -> float:
        """Calcula o percentual de bases Guanina (G) e Citosina (C)."""
        total_length = len(self)
        if total_length == 0:
            return 0.0
        
        comp = self.get_base_composition()
        gc_count = comp.get('G', 0) + comp.get('C', 0)
        return (gc_count / total_length) * 100.0
    
    # --- COMPLEMENTARIDADE ---
    
    def complement(self) -> 'DNASequence':
        """Retorna a sequência complementar."""
        complement_seq = self._sequence.translate(self._COMPLEMENT_MAP)
        return DNASequence(complement_seq, f"{self.id}_complement")
    
    def reverse(self) -> 'DNASequence':
        """Retorna a sequência invertida (não o complemento)."""
        return DNASequence(
            self._sequence[::-1],
            f"{self.id}_reverse"
        )
    
    def reverse_complement(self) -> 'DNASequence':
        """Retorna o complemento reverso da sequência."""
        rev_comp_seq = self.complement().sequence[::-1]
        return DNASequence(rev_comp_seq, f"{self.id}_revcomp")
    
    # --- TRANSCRIÇÃO ---
    
    def transcribe(self) -> str:
        """Transcreve a sequência de DNA para RNA (substituindo T por U)."""
        return self._sequence.replace('T', 'U')
    
    # --- BUSCA DE PADRÕES ---
    
    def find_all_occurrences(self, pattern: str) -> List[int]:
        """Encontra todas as posições iniciais (índices) de um padrão."""
        target_pattern = pattern.upper()
        positions = []
        start_index = 0
        
        while True:
            position = self._sequence.find(target_pattern, start_index)
            if position == -1:
                break
            
            positions.append(position)
            start_index = position + 1
            
        return positions
    
    def count_pattern_occurrences(self, pattern: str) -> int:
        """Conta o número de ocorrências exatas de um padrão."""
        return self._sequence.count(pattern.upper())
    
    # --- ANÁLISE DE ORFs ---
    
    def find_orfs(self, min_length_bp: int = 100) -> List[Tuple[int, int, str]]:
        """
        Identifica possíveis Open Reading Frames (ORFs) na direção forward.
        Retorna uma lista de tuplas: (start_index, end_index_exclusive, orf_sequence).
        """
        orfs = []
        START_CODON = 'ATG'
        STOP_CODONS = {'TAA', 'TAG', 'TGA'}
        
        for frame_shift in range(3):
            i = frame_shift
            while i < len(self) - 2:
                codon = self._sequence[i:i+3]
                
                if codon == START_CODON:
                    start_pos = i
                    j = i + 3
                    
                    while j < len(self) - 2:
                        stop_codon = self._sequence[j:j+3]
                        
                        if stop_codon in STOP_CODONS:
                            end_pos_exclusive = j + 3
                            orf_length = end_pos_exclusive - start_pos
                            
                            if orf_length >= min_length_bp:
                                orf_seq = self._sequence[start_pos:end_pos_exclusive]
                                orfs.append((start_pos, end_pos_exclusive, orf_seq))
                            
                            i = end_pos_exclusive
                            break
                        
                        j += 3
                    else:
                        i += 3
                else:
                    i += 3
                    
        return orfs

    # --- ANÁLISES AVANÇADAS ---
    
    def calculate_melting_temp(self) -> float:
        """
        Calcula a temperatura de desnaturação (Tm) da sequência.
        """
        comp = self.get_base_composition()
        gc_count = comp.get('G', 0) + comp.get('C', 0)
        at_count = comp.get('A', 0) + comp.get('T', 0)
        
        total_length = len(self)
        if total_length == 0:
            return 0.0

        if total_length < 14:
            # Regra de Wallace (rápida)
            return 4.0 * gc_count + 2.0 * at_count
        else:
            # Fórmula empírica para sequências mais longas
            gc_fraction = self.gc_content() / 100.0
            return 64.9 + 41.0 * (gc_fraction - (16.4 / total_length))
    
    def find_restriction_sites(self, enzyme: str) -> List[int]:
        """
        Encontra as posições dos sítios de restrição de uma enzima comum.
        """
        RESTRICTION_SITES = {
            'EcoRI': 'GAATTC',
            'BamHI': 'GGATCC',
            'HindIII': 'AAGCTT',
            'PstI': 'CTGCAG',
            'SmaI': 'CCCGGG',
            'XbaI': 'TCTAGA'
        }
        
        recognition_sequence = RESTRICTION_SITES.get(enzyme)
        
        if recognition_sequence is None:
            return []
        
        return self.find_all_occurrences(recognition_sequence)


if __name__ == "__main__":
    dna_seq = DNASequence("ATGCATGC", seq_id="exemplo_teste")
    print(f"ID: {dna_seq.id}")
    print(f"GC%: {dna_seq.gc_content():.2f}%")