package util

import (
	"fmt"
	"strings"
)

// 密码子表 (标准遗传密码)
var CodonTable = map[string]string{
	"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
	"CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
	"ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
	"GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
	"TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
	"CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
	"ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
	"GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
	"TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
	"CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
	"AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
	"GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
	"TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
	"CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
	"AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
	"GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

// DNA转氨基酸序列
func DnaToProtein(dna string) (string, error) {
	dna = strings.ToUpper(strings.ReplaceAll(dna, " ", ""))

	// 检查DNA长度是否为3的倍数
	if len(dna)%3 != 0 {
		return "", fmt.Errorf("DNA序列长度不是3的倍数")
	}

	var protein strings.Builder

	for i := 0; i < len(dna); i += 3 {
		codon := dna[i : i+3]
		aa, ok := CodonTable[codon]
		if !ok {
			return "", fmt.Errorf("无效的密码子: %s", codon)
		}
		protein.WriteString(aa)
	}

	return protein.String(), nil
}
