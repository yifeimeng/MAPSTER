V29 :0x4 cufft_wrapper
21 mod_CUFFT_wrapper.f90 S624 0
03/07/2017  12:04:36
use cufft public 0 direct
use iso_c_binding private
enduse
D 56 24 644 8 643 7
D 62 24 646 8 645 7
D 130 21 12 1 3 113 0 0 1 0 0
 0 112 3 3 113 113
D 133 21 12 1 3 113 0 0 1 0 0
 0 112 3 3 113 113
D 136 21 12 1 3 115 0 0 1 0 0
 0 114 3 3 115 115
D 139 21 12 1 3 115 0 0 1 0 0
 0 114 3 3 115 115
D 142 21 12 2 116 121 0 0 1 0 0
 0 117 3 3 118 118
 0 119 118 3 120 120
D 145 21 12 2 116 121 0 0 1 0 0
 0 117 3 3 118 118
 0 119 118 3 120 120
D 148 21 12 2 122 127 0 0 1 0 0
 0 123 3 3 124 124
 0 125 124 3 126 126
D 151 21 12 2 122 127 0 0 1 0 0
 0 123 3 3 124 124
 0 125 124 3 126 126
D 154 21 12 3 128 136 0 0 1 0 0
 0 129 3 3 130 130
 0 131 130 3 132 132
 0 133 134 3 135 135
D 157 21 12 3 128 136 0 0 1 0 0
 0 129 3 3 130 130
 0 131 130 3 132 132
 0 133 134 3 135 135
D 160 21 12 3 137 145 0 0 1 0 0
 0 138 3 3 139 139
 0 140 139 3 141 141
 0 142 143 3 144 144
D 163 21 12 3 137 145 0 0 1 0 0
 0 138 3 3 139 139
 0 140 139 3 141 141
 0 142 143 3 144 144
D 166 21 11 1 3 147 0 0 1 0 0
 0 146 3 3 147 147
D 169 21 11 1 3 147 0 0 1 0 0
 0 146 3 3 147 147
D 172 21 11 1 3 149 0 0 1 0 0
 0 148 3 3 149 149
D 175 21 11 1 3 149 0 0 1 0 0
 0 148 3 3 149 149
D 178 21 11 2 150 155 0 0 1 0 0
 0 151 3 3 152 152
 0 153 152 3 154 154
D 181 21 11 2 150 155 0 0 1 0 0
 0 151 3 3 152 152
 0 153 152 3 154 154
D 184 21 11 2 156 161 0 0 1 0 0
 0 157 3 3 158 158
 0 159 158 3 160 160
D 187 21 11 2 156 161 0 0 1 0 0
 0 157 3 3 158 158
 0 159 158 3 160 160
D 190 21 11 3 162 170 0 0 1 0 0
 0 163 3 3 164 164
 0 165 164 3 166 166
 0 167 168 3 169 169
D 193 21 11 3 162 170 0 0 1 0 0
 0 163 3 3 164 164
 0 165 164 3 166 166
 0 167 168 3 169 169
D 196 21 11 3 171 179 0 0 1 0 0
 0 172 3 3 173 173
 0 174 173 3 175 175
 0 176 177 3 178 178
D 199 21 11 3 171 179 0 0 1 0 0
 0 172 3 3 173 173
 0 174 173 3 175 175
 0 176 177 3 178 178
S 624 24 0 0 0 8 1 0 5011 10005 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20 0 0 0 0 0 0 cufft_wrapper
R 643 25 6 iso_c_binding c_ptr
R 644 5 7 iso_c_binding val c_ptr
R 645 25 8 iso_c_binding c_funptr
R 646 5 9 iso_c_binding val c_funptr
R 679 6 42 iso_c_binding c_null_ptr$ac
R 681 6 44 iso_c_binding c_null_funptr$ac
R 682 26 45 iso_c_binding ==
R 684 26 47 iso_c_binding !=
S 811 19 0 0 0 8 1 624 6309 4000 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 22 2 0 0 0 0 0 624 0 0 0 0 fft1
O 811 2 813 812
S 812 27 0 0 0 8 829 624 6314 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 43 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 dfft1d
Q 812 811 0
S 813 27 0 0 0 8 903 624 6321 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 49 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 sfft1d
Q 813 811 0
S 814 19 0 0 0 6 1 624 6328 4000 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 26 2 0 0 0 0 0 624 0 0 0 0 ifft1
O 814 2 816 815
S 815 27 0 0 0 8 837 624 6334 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 44 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 dfft1b
Q 815 814 0
S 816 27 0 0 0 8 911 624 6341 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 50 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 sfft1b
Q 816 814 0
S 817 19 0 0 0 8 1 624 6348 4000 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 30 2 0 0 0 0 0 624 0 0 0 0 fft2
O 817 2 819 818
S 818 27 0 0 0 8 845 624 6353 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 45 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 dfft2d
Q 818 817 0
S 819 27 0 0 0 8 919 624 6360 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 51 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 sfft2d
Q 819 817 0
S 820 19 0 0 0 6 1 624 6367 4000 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 34 2 0 0 0 0 0 624 0 0 0 0 ifft2
O 820 2 822 821
S 821 27 0 0 0 8 857 624 6373 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 46 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 dfft2b
Q 821 820 0
S 822 27 0 0 0 8 931 624 6380 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 52 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 sfft2b
Q 822 820 0
S 823 19 0 0 0 8 1 624 6387 4000 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 38 2 0 0 0 0 0 624 0 0 0 0 fft3
O 823 2 825 824
S 824 27 0 0 0 8 869 624 6392 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 47 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 dfft3d
Q 824 823 0
S 825 27 0 0 0 8 943 624 6399 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 53 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 sfft3d
Q 825 823 0
S 826 19 0 0 0 6 1 624 6406 4000 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 42 2 0 0 0 0 0 624 0 0 0 0 ifft3
O 826 2 828 827
S 827 27 0 0 0 8 886 624 6412 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 48 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 dfft3b
Q 827 826 0
S 828 27 0 0 0 8 960 624 6419 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 54 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 sfft3b
Q 828 826 0
S 829 23 5 0 0 0 835 624 6314 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 dfft1d
S 830 6 3 0 0 6 1 829 6426 800004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiyb
S 831 7 3 0 0 130 1 829 6433 800204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array_in
S 832 1 3 0 0 6 1 829 6442 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiy
S 833 7 3 0 0 133 1 829 6448 800204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array_out
S 834 1 3 0 0 6 1 829 6458 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopix
S 835 14 5 0 0 0 1 829 6314 200 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 61 5 0 0 0 0 0 0 0 0 0 0 0 0 59 0 624 0 0 0 0 dfft1d
F 835 5 830 831 832 833 834
S 836 6 1 0 0 6 1 829 6464 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_112
S 837 23 5 0 0 0 843 624 6334 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 dfft1b
S 838 6 3 0 0 6 1 837 6426 800004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiyb
S 839 7 3 0 0 136 1 837 6433 800204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array_in
S 840 1 3 0 0 6 1 837 6442 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiy
S 841 7 3 0 0 139 1 837 6448 800204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array_out
S 842 1 3 0 0 6 1 837 6458 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopix
S 843 14 5 0 0 0 1 837 6334 200 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 67 5 0 0 0 0 0 0 0 0 0 0 0 0 93 0 624 0 0 0 0 dfft1b
F 843 5 838 839 840 841 842
S 844 6 1 0 0 6 1 837 6472 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_114
S 845 23 5 0 0 0 852 624 6353 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 dfft2d
S 846 6 3 0 0 6 1 845 6426 800004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiyb
S 847 6 3 0 0 6 1 845 6480 800004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopixb
S 848 7 3 0 0 142 1 845 6433 800204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array_in
S 849 1 3 0 0 6 1 845 6442 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiy
S 850 7 3 0 0 145 1 845 6448 800204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array_out
S 851 1 3 0 0 6 1 845 6458 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopix
S 852 14 5 0 0 0 1 845 6353 200 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 73 6 0 0 0 0 0 0 0 0 0 0 0 0 133 0 624 0 0 0 0 dfft2d
F 852 6 846 847 848 849 850 851
S 853 6 1 0 0 6 1 845 6487 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_116
S 854 6 1 0 0 6 1 845 6495 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_117
S 855 6 1 0 0 6 1 845 6503 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_120
S 856 6 1 0 0 6 1 845 6511 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_122
S 857 23 5 0 0 0 864 624 6373 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 dfft2b
S 858 6 3 0 0 6 1 857 6426 800004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiyb
S 859 6 3 0 0 6 1 857 6480 800004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopixb
S 860 7 3 0 0 148 1 857 6433 800204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array_in
S 861 1 3 0 0 6 1 857 6442 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiy
S 862 7 3 0 0 151 1 857 6448 800204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array_out
S 863 1 3 0 0 6 1 857 6458 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopix
S 864 14 5 0 0 0 1 857 6373 200 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 80 6 0 0 0 0 0 0 0 0 0 0 0 0 167 0 624 0 0 0 0 dfft2b
F 864 6 858 859 860 861 862 863
S 865 6 1 0 0 6 1 857 6511 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_122
S 866 6 1 0 0 6 1 857 6519 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_123
S 867 6 1 0 0 6 1 857 6527 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_126
S 868 6 1 0 0 6 1 857 6535 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_128
S 869 23 5 0 0 0 879 624 6392 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 dfft3d
S 870 6 3 0 0 6 1 869 6426 800004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiyb
S 871 6 3 0 0 6 1 869 6480 800004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopixb
S 872 6 3 0 0 6 1 869 6543 800004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopizb
S 873 7 3 0 0 154 1 869 6433 800204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array_in
S 874 1 3 0 0 6 1 869 6442 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiy
S 875 1 3 0 0 6 1 869 6458 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopix
S 876 7 3 0 0 157 1 869 6448 800204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array_out
S 877 1 3 0 0 6 1 869 6550 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiya
S 878 1 3 0 0 6 1 869 6557 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopixa
S 879 14 5 0 0 0 1 869 6392 200 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 87 9 0 0 0 0 0 0 0 0 0 0 0 0 201 0 624 0 0 0 0 dfft3d
F 879 9 870 871 872 873 874 875 876 877 878
S 880 6 1 0 0 6 1 869 6535 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_128
S 881 6 1 0 0 6 1 869 6564 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_129
S 882 6 1 0 0 6 1 869 6572 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_133
S 883 6 1 0 0 6 1 869 6580 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_130
S 884 6 1 0 0 6 1 869 6588 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_136
S 885 6 1 0 0 6 1 869 6596 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_139
S 886 23 5 0 0 0 896 624 6412 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 dfft3b
S 887 6 3 0 0 6 1 886 6426 800004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiyb
S 888 6 3 0 0 6 1 886 6480 800004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopixb
S 889 6 3 0 0 6 1 886 6543 800004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopizb
S 890 7 3 0 0 160 1 886 6433 800204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array_in
S 891 1 3 0 0 6 1 886 6442 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiy
S 892 1 3 0 0 6 1 886 6458 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopix
S 893 7 3 0 0 163 1 886 6448 800204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array_out
S 894 1 3 0 0 6 1 886 6550 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiya
S 895 1 3 0 0 6 1 886 6557 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopixa
S 896 14 5 0 0 0 1 886 6412 200 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 97 9 0 0 0 0 0 0 0 0 0 0 0 0 235 0 624 0 0 0 0 dfft3b
F 896 9 887 888 889 890 891 892 893 894 895
S 897 6 1 0 0 6 1 886 6604 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_137
S 898 6 1 0 0 6 1 886 6612 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_138
S 899 6 1 0 0 6 1 886 6620 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_142
S 900 6 1 0 0 6 1 886 6596 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_139
S 901 6 1 0 0 6 1 886 6628 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_145
S 902 6 1 0 0 6 1 886 6636 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_148
S 903 23 5 0 0 0 909 624 6321 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 sfft1d
S 904 6 3 0 0 6 1 903 6426 800004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiyb
S 905 7 3 0 0 166 1 903 6433 800204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array_in
S 906 1 3 0 0 6 1 903 6442 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiy
S 907 7 3 0 0 169 1 903 6448 800204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array_out
S 908 1 3 0 0 6 1 903 6458 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopix
S 909 14 5 0 0 0 1 903 6321 200 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 107 5 0 0 0 0 0 0 0 0 0 0 0 0 273 0 624 0 0 0 0 sfft1d
F 909 5 904 905 906 907 908
S 910 6 1 0 0 6 1 903 6644 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_146
S 911 23 5 0 0 0 917 624 6341 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 sfft1b
S 912 6 3 0 0 6 1 911 6426 800004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiyb
S 913 7 3 0 0 172 1 911 6433 800204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array_in
S 914 1 3 0 0 6 1 911 6442 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiy
S 915 7 3 0 0 175 1 911 6448 800204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array_out
S 916 1 3 0 0 6 1 911 6458 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopix
S 917 14 5 0 0 0 1 911 6341 200 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 113 5 0 0 0 0 0 0 0 0 0 0 0 0 306 0 624 0 0 0 0 sfft1b
F 917 5 912 913 914 915 916
S 918 6 1 0 0 6 1 911 6636 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_148
S 919 23 5 0 0 0 926 624 6360 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 sfft2d
S 920 6 3 0 0 6 1 919 6426 800004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiyb
S 921 6 3 0 0 6 1 919 6480 800004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopixb
S 922 7 3 0 0 178 1 919 6433 800204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array_in
S 923 1 3 0 0 6 1 919 6442 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiy
S 924 7 3 0 0 181 1 919 6448 800204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array_out
S 925 1 3 0 0 6 1 919 6458 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopix
S 926 14 5 0 0 0 1 919 6360 200 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 119 6 0 0 0 0 0 0 0 0 0 0 0 0 344 0 624 0 0 0 0 sfft2d
F 926 6 920 921 922 923 924 925
S 927 6 1 0 0 6 1 919 6652 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_150
S 928 6 1 0 0 6 1 919 6660 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_151
S 929 6 1 0 0 6 1 919 6668 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_154
S 930 6 1 0 0 6 1 919 6676 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_156
S 931 23 5 0 0 0 938 624 6380 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 sfft2b
S 932 6 3 0 0 6 1 931 6426 800004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiyb
S 933 6 3 0 0 6 1 931 6480 800004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopixb
S 934 7 3 0 0 184 1 931 6433 800204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array_in
S 935 1 3 0 0 6 1 931 6442 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiy
S 936 7 3 0 0 187 1 931 6448 800204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array_out
S 937 1 3 0 0 6 1 931 6458 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopix
S 938 14 5 0 0 0 1 931 6380 200 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 126 6 0 0 0 0 0 0 0 0 0 0 0 0 377 0 624 0 0 0 0 sfft2b
F 938 6 932 933 934 935 936 937
S 939 6 1 0 0 6 1 931 6676 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_156
S 940 6 1 0 0 6 1 931 6684 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_157
S 941 6 1 0 0 6 1 931 6692 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_160
S 942 6 1 0 0 6 1 931 6700 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_162
S 943 23 5 0 0 0 953 624 6399 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 sfft3d
S 944 6 3 0 0 6 1 943 6426 800004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiyb
S 945 6 3 0 0 6 1 943 6480 800004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopixb
S 946 6 3 0 0 6 1 943 6543 800004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopizb
S 947 7 3 0 0 190 1 943 6433 800204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array_in
S 948 1 3 0 0 6 1 943 6442 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiy
S 949 1 3 0 0 6 1 943 6458 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopix
S 950 7 3 0 0 193 1 943 6448 800204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array_out
S 951 1 3 0 0 6 1 943 6550 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiya
S 952 1 3 0 0 6 1 943 6557 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopixa
S 953 14 5 0 0 0 1 943 6399 200 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 133 9 0 0 0 0 0 0 0 0 0 0 0 0 410 0 624 0 0 0 0 sfft3d
F 953 9 944 945 946 947 948 949 950 951 952
S 954 6 1 0 0 6 1 943 6700 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_162
S 955 6 1 0 0 6 1 943 6708 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_163
S 956 6 1 0 0 6 1 943 6716 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_167
S 957 6 1 0 0 6 1 943 6724 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_164
S 958 6 1 0 0 6 1 943 6732 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_170
S 959 6 1 0 0 6 1 943 6740 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_173
S 960 23 5 0 0 0 970 624 6419 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 sfft3b
S 961 6 3 0 0 6 1 960 6426 800004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiyb
S 962 6 3 0 0 6 1 960 6480 800004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopixb
S 963 6 3 0 0 6 1 960 6543 800004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopizb
S 964 7 3 0 0 196 1 960 6433 800204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array_in
S 965 1 3 0 0 6 1 960 6442 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiy
S 966 1 3 0 0 6 1 960 6458 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopix
S 967 7 3 0 0 199 1 960 6448 800204 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array_out
S 968 1 3 0 0 6 1 960 6550 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopiya
S 969 1 3 0 0 6 1 960 6557 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 nopixa
S 970 14 5 0 0 0 1 960 6419 200 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 143 9 0 0 0 0 0 0 0 0 0 0 0 0 444 0 624 0 0 0 0 sfft3b
F 970 9 961 962 963 964 965 966 967 968 969
S 971 6 1 0 0 6 1 960 6748 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_171
S 972 6 1 0 0 6 1 960 6756 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_172
S 973 6 1 0 0 6 1 960 6764 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_176
S 974 6 1 0 0 6 1 960 6740 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_173
S 975 6 1 0 0 6 1 960 6772 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_179
S 976 6 1 0 0 6 1 960 6780 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_182
A 67 1 0 0 0 56 679 0 0 0 0 0 0 0 0 0 0 0 0 0
A 70 1 0 0 0 62 681 0 0 0 0 0 0 0 0 0 0 0 0 0
A 112 1 0 0 108 6 830 0 0 0 0 0 0 0 0 0 0 0 0 0
A 113 1 0 0 72 6 836 0 0 0 0 0 0 0 0 0 0 0 0 0
A 114 1 0 0 0 6 838 0 0 0 0 0 0 0 0 0 0 0 0 0
A 115 1 0 0 0 6 844 0 0 0 0 0 0 0 0 0 0 0 0 0
A 116 1 0 0 48 6 856 0 0 0 0 0 0 0 0 0 0 0 0 0
A 117 1 0 0 0 6 846 0 0 0 0 0 0 0 0 0 0 0 0 0
A 118 1 0 0 0 6 853 0 0 0 0 0 0 0 0 0 0 0 0 0
A 119 1 0 0 9 6 847 0 0 0 0 0 0 0 0 0 0 0 0 0
A 120 1 0 0 44 6 854 0 0 0 0 0 0 0 0 0 0 0 0 0
A 121 1 0 0 46 6 855 0 0 0 0 0 0 0 0 0 0 0 0 0
A 122 1 0 0 0 6 868 0 0 0 0 0 0 0 0 0 0 0 0 0
A 123 1 0 0 52 6 858 0 0 0 0 0 0 0 0 0 0 0 0 0
A 124 1 0 0 0 6 865 0 0 0 0 0 0 0 0 0 0 0 0 0
A 125 1 0 0 54 6 859 0 0 0 0 0 0 0 0 0 0 0 0 0
A 126 1 0 0 0 6 866 0 0 0 0 0 0 0 0 0 0 0 0 0
A 127 1 0 0 0 6 867 0 0 0 0 0 0 0 0 0 0 0 0 0
A 128 1 0 0 0 6 885 0 0 0 0 0 0 0 0 0 0 0 0 0
A 129 1 0 0 0 6 870 0 0 0 0 0 0 0 0 0 0 0 0 0
A 130 1 0 0 0 6 880 0 0 0 0 0 0 0 0 0 0 0 0 0
A 131 1 0 0 0 6 871 0 0 0 0 0 0 0 0 0 0 0 0 0
A 132 1 0 0 0 6 881 0 0 0 0 0 0 0 0 0 0 0 0 0
A 133 1 0 0 0 6 872 0 0 0 0 0 0 0 0 0 0 0 0 0
A 134 1 0 0 0 6 882 0 0 0 0 0 0 0 0 0 0 0 0 0
A 135 1 0 0 0 6 883 0 0 0 0 0 0 0 0 0 0 0 0 0
A 136 1 0 0 0 6 884 0 0 0 0 0 0 0 0 0 0 0 0 0
A 137 1 0 0 0 6 902 0 0 0 0 0 0 0 0 0 0 0 0 0
A 138 1 0 0 0 6 887 0 0 0 0 0 0 0 0 0 0 0 0 0
A 139 1 0 0 0 6 897 0 0 0 0 0 0 0 0 0 0 0 0 0
A 140 1 0 0 0 6 888 0 0 0 0 0 0 0 0 0 0 0 0 0
A 141 1 0 0 0 6 898 0 0 0 0 0 0 0 0 0 0 0 0 0
A 142 1 0 0 0 6 889 0 0 0 0 0 0 0 0 0 0 0 0 0
A 143 1 0 0 0 6 899 0 0 0 0 0 0 0 0 0 0 0 0 0
A 144 1 0 0 0 6 900 0 0 0 0 0 0 0 0 0 0 0 0 0
A 145 1 0 0 0 6 901 0 0 0 0 0 0 0 0 0 0 0 0 0
A 146 1 0 0 0 6 904 0 0 0 0 0 0 0 0 0 0 0 0 0
A 147 1 0 0 0 6 910 0 0 0 0 0 0 0 0 0 0 0 0 0
A 148 1 0 0 0 6 912 0 0 0 0 0 0 0 0 0 0 0 0 0
A 149 1 0 0 0 6 918 0 0 0 0 0 0 0 0 0 0 0 0 0
A 150 1 0 0 0 6 930 0 0 0 0 0 0 0 0 0 0 0 0 0
A 151 1 0 0 0 6 920 0 0 0 0 0 0 0 0 0 0 0 0 0
A 152 1 0 0 0 6 927 0 0 0 0 0 0 0 0 0 0 0 0 0
A 153 1 0 0 0 6 921 0 0 0 0 0 0 0 0 0 0 0 0 0
A 154 1 0 0 0 6 928 0 0 0 0 0 0 0 0 0 0 0 0 0
A 155 1 0 0 0 6 929 0 0 0 0 0 0 0 0 0 0 0 0 0
A 156 1 0 0 0 6 942 0 0 0 0 0 0 0 0 0 0 0 0 0
A 157 1 0 0 0 6 932 0 0 0 0 0 0 0 0 0 0 0 0 0
A 158 1 0 0 0 6 939 0 0 0 0 0 0 0 0 0 0 0 0 0
A 159 1 0 0 0 6 933 0 0 0 0 0 0 0 0 0 0 0 0 0
A 160 1 0 0 0 6 940 0 0 0 0 0 0 0 0 0 0 0 0 0
A 161 1 0 0 0 6 941 0 0 0 0 0 0 0 0 0 0 0 0 0
A 162 1 0 0 0 6 959 0 0 0 0 0 0 0 0 0 0 0 0 0
A 163 1 0 0 0 6 944 0 0 0 0 0 0 0 0 0 0 0 0 0
A 164 1 0 0 0 6 954 0 0 0 0 0 0 0 0 0 0 0 0 0
A 165 1 0 0 0 6 945 0 0 0 0 0 0 0 0 0 0 0 0 0
A 166 1 0 0 5 6 955 0 0 0 0 0 0 0 0 0 0 0 0 0
A 167 1 0 0 0 6 946 0 0 0 0 0 0 0 0 0 0 0 0 0
A 168 1 0 0 0 6 956 0 0 0 0 0 0 0 0 0 0 0 0 0
A 169 1 0 0 7 6 957 0 0 0 0 0 0 0 0 0 0 0 0 0
A 170 1 0 0 0 6 958 0 0 0 0 0 0 0 0 0 0 0 0 0
A 171 1 0 0 0 6 976 0 0 0 0 0 0 0 0 0 0 0 0 0
A 172 1 0 0 0 6 961 0 0 0 0 0 0 0 0 0 0 0 0 0
A 173 1 0 0 0 6 971 0 0 0 0 0 0 0 0 0 0 0 0 0
A 174 1 0 0 0 6 962 0 0 0 0 0 0 0 0 0 0 0 0 0
A 175 1 0 0 0 6 972 0 0 0 0 0 0 0 0 0 0 0 0 0
A 176 1 0 0 0 6 963 0 0 0 0 0 0 0 0 0 0 0 0 0
A 177 1 0 0 0 6 973 0 0 0 0 0 0 0 0 0 0 0 0 0
A 178 1 0 0 0 6 974 0 0 0 0 0 0 0 0 0 0 0 0 0
A 179 1 0 0 0 6 975 0 0 0 0 0 0 0 0 0 0 0 0 0
Z
J 149 1 1
V 67 56 7 0
S 0 56 0 0 0
A 0 6 0 0 1 2 0
J 150 1 1
V 70 62 7 0
S 0 62 0 0 0
A 0 6 0 0 1 2 0
Z