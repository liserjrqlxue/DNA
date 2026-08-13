package util

import (
	"bytes"
	"fmt"
	"testing"
)

// ---------- 单元测试：验证所有实现正确 ----------
func TestReverseComplementFunctions(t *testing.T) {
	tests := []struct {
		input    string
		expected string
	}{
		{"", ""},
		{"A", "T"},
		{"T", "A"},
		{"C", "G"},
		{"G", "C"},
		{"a", "t"},
		{"t", "a"},
		{"c", "g"},
		{"g", "c"},
		{"ATCG", "CGAT"},
		{"atcg", "cgat"},
		{"AaTtCcGg", "cCgGaAtT"},
		{"ATCGN", "NCGAT"},
		{"XYZ", "ZYX"},
		{"ACGTacgtNn", "nNacgtACGT"},
	}

	for _, tt := range tests {
		// 测试 ReverseComplement (string -> string)
		if got := ReverseComplement(tt.input); got != tt.expected {
			t.Errorf("ReverseComplement(%q) = %q, want %q", tt.input, got, tt.expected)
		}

		// 测试 ReverseComplement2Bytes (string -> []byte)
		gotBytes := ReverseComplement2Bytes(tt.input)
		if !bytes.Equal(gotBytes, []byte(tt.expected)) {
			t.Errorf("ReverseComplement2Bytes(%q) = %q, want %q", tt.input, gotBytes, tt.expected)
		}

		// 测试 ReverseComplementBytes ([]byte -> []byte)
		inputBytes := []byte(tt.input)
		gotBytes2 := ReverseComplementBytes(inputBytes)
		if !bytes.Equal(gotBytes2, []byte(tt.expected)) {
			t.Errorf("ReverseComplementBytes(%q) = %q, want %q", tt.input, gotBytes2, tt.expected)
		}
	}
}

// ---------- 基准测试：按输入类型和长度分组 ----------
var seqSizes = []int{10, 100, 1000, 10000, 100000}

// 辅助函数：生成指定长度的 DNA 序列（带 N 和大小写混合）
func makeTestSeq(size int) (string, []byte) {
	// 用 bytes.Repeat 生成 “ATCGN” 的循环，以包含 N 和大小写
	base := []byte("ATCGatcgNn")
	repeat := (size + len(base) - 1) / len(base)
	b := bytes.Repeat(base, repeat)[:size]
	return string(b), b
}

// 辅助函数：itoa（避免 import strconv）
func itoa(i int) string {
	if i == 0 {
		return "0"
	}
	var buf [20]byte
	pos := len(buf)
	for i > 0 {
		pos--
		buf[pos] = byte('0' + i%10)
		i /= 10
	}
	return string(buf[pos:])
}

// BenchmarkReverseComplement 按长度和实现方式分组
func BenchmarkReverseComplement(b *testing.B) {
	for _, size := range seqSizes {
		s, seq := makeTestSeq(size)
		b.Run("size="+itoa(size), func(b *testing.B) {
			// ---- 输入为 []byte ----
			b.Run("BytesInput/ReverseComplementBytes→[]byte", func(b *testing.B) {
				b.SetBytes(int64(size))
				for i := 0; i < b.N; i++ {
					_ = ReverseComplementBytes(seq)
				}
			})
			b.Run("BytesInput/ReverseComplementBytes→string", func(b *testing.B) {
				b.SetBytes(int64(size))
				for i := 0; i < b.N; i++ {
					_ = string(ReverseComplementBytes(seq))
				}
			})

			// ---- 输入为 string ----
			b.Run("StringInput/ReverseComplement2Bytes→[]byte", func(b *testing.B) {
				b.SetBytes(int64(size))
				for i := 0; i < b.N; i++ {
					_ = ReverseComplement2Bytes(s)
				}
			})
			b.Run("StringInput/ReverseComplement2Bytes→string", func(b *testing.B) {
				b.SetBytes(int64(size))
				for i := 0; i < b.N; i++ {
					_ = string(ReverseComplement2Bytes(s))
				}
			})
			b.Run("StringInput/ReverseComplement→string", func(b *testing.B) {
				b.SetBytes(int64(size))
				for i := 0; i < b.N; i++ {
					_ = ReverseComplement(s)
				}
			})
		})
	}
}

// ---------- 示例测试：演示用法（使用 fmt.Println 写入 stdout） ----------
func ExampleReverseComplement() {
	s := "ATCGN"
	result := ReverseComplement(s)
	fmt.Println(result)
	// Output: NCGAT
}

func ExampleReverseComplementBytes() {
	seq := []byte("ATCGN")
	result := ReverseComplementBytes(seq)
	fmt.Println(string(result))
	// Output: NCGAT
}

func ExampleReverseComplement2Bytes() {
	s := "ATCGN"
	result := ReverseComplement2Bytes(s)
	fmt.Println(string(result))
	// Output: NCGAT
}
