V29 :0x4 m_tilt
10 m_tilt.f90 S624 0
03/07/2017  13:36:38
use global_variables private
use m_precision private
enduse
D 180 21 12 2 231 237 1 1 0 0 1
 3 232 3 3 232 233
 3 234 235 3 234 236
D 183 21 12 2 238 244 1 1 0 0 1
 3 239 3 3 239 240
 3 241 242 3 241 243
D 186 21 12 2 245 251 1 1 0 0 1
 3 246 3 3 246 247
 3 248 249 3 248 250
D 189 21 6 1 3 253 0 0 1 0 0
 0 252 3 3 253 253
D 192 21 9 2 254 260 1 1 0 0 1
 3 255 3 3 255 256
 3 257 258 3 257 259
D 195 21 9 2 261 269 0 0 1 0 0
 0 264 3 3 265 265
 0 267 265 3 268 268
S 624 24 0 0 0 6 1 0 5011 10005 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 m_tilt
S 626 23 0 0 0 8 634 624 5030 4 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 fp_kind
S 628 23 0 0 0 6 671 624 5055 4 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 ifactory
S 629 23 0 0 0 6 670 624 5064 4 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 ifactorx
R 634 16 3 m_precision fp_kind
R 670 6 10 global_variables ifactorx
R 671 6 11 global_variables ifactory
S 822 6 4 0 0 16 1 624 6383 80000c 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 826 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 tilt_illumination
S 824 6 4 0 0 9 825 624 6401 4 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 827 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 alpha
S 825 6 4 0 0 9 1 624 6407 4 0 A 0 0 0 0 B 0 0 0 0 0 8 0 0 0 0 0 0 827 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 beta
S 826 11 0 0 0 8 819 624 6412 40800000 805000 A 0 0 0 0 B 0 0 0 0 0 4 0 0 822 822 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _m_tilt$8
S 827 11 0 0 0 8 826 624 6422 40800000 805000 A 0 0 0 0 B 0 0 0 0 0 16 0 0 824 825 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _m_tilt$2
S 828 23 5 0 0 0 829 624 6432 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 prompt_tilt
S 829 14 5 0 0 0 1 828 6432 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 14 0 624 0 0 0 0 prompt_tilt
F 829 0
S 830 23 5 0 0 0 831 624 6444 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 setup_tilt
S 831 14 5 0 0 0 1 830 6444 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 53 0 624 0 0 0 0 setup_tilt
F 831 0
S 832 23 5 0 0 0 833 624 6455 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 setup_specimen_tilt
S 833 14 5 0 0 0 1 832 6455 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 5 0 0 0 0 0 0 0 0 0 0 0 0 0 96 0 624 0 0 0 0 setup_specimen_tilt
F 833 0
S 834 23 5 0 0 0 836 624 6475 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 tilt_wave_function
S 835 7 3 0 0 180 1 834 6494 20000004 10003000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 psi
S 836 14 5 0 0 0 1 834 6475 20000000 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 6 1 0 0 0 0 0 0 0 0 0 0 0 0 125 0 624 0 0 0 0 tilt_wave_function
F 836 1 835
S 837 6 1 0 0 6 1 834 6498 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_0_1
S 838 6 1 0 0 6 1 834 6506 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_2_1
S 839 6 1 0 0 6 1 834 6514 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_3_1
S 840 6 1 0 0 6 1 834 6522 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_5_1
S 841 6 1 0 0 6 1 834 6530 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_6_1
S 842 6 1 0 0 6 1 834 6538 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_239
S 843 6 1 0 0 6 1 834 6546 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_242
S 844 23 5 0 0 0 846 624 6554 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 tilt_wave_function_phase_ramp
S 845 7 3 0 0 183 1 844 6494 20000004 10003000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 psi
S 846 14 5 0 0 0 1 844 6554 20000000 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 8 1 0 0 0 0 0 0 0 0 0 0 0 0 138 0 624 0 0 0 0 tilt_wave_function_phase_ramp
F 846 1 845
S 847 6 1 0 0 6 1 844 6498 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_0_1
S 848 6 1 0 0 6 1 844 6506 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_2_1
S 849 6 1 0 0 6 1 844 6514 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_3_1
S 850 6 1 0 0 6 1 844 6522 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_5_1
S 851 6 1 0 0 6 1 844 6530 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_6_1
S 852 6 1 0 0 6 1 844 6584 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_310
S 853 6 1 0 0 6 1 844 6592 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_313
S 854 23 5 0 0 0 856 624 6600 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 tilt_wave_function_shift
S 855 7 3 0 0 186 1 854 6494 20000004 10003000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 psi
S 856 14 5 0 0 0 1 854 6600 20000000 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 10 1 0 0 0 0 0 0 0 0 0 0 0 0 170 0 624 0 0 0 0 tilt_wave_function_shift
F 856 1 855
S 857 6 1 0 0 6 1 854 6498 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_0_1
S 858 6 1 0 0 6 1 854 6506 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_2_1
S 859 6 1 0 0 6 1 854 6514 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_3_1
S 860 6 1 0 0 6 1 854 6522 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_5_1
S 861 6 1 0 0 6 1 854 6530 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_6_1
S 862 6 1 0 0 6 1 854 6625 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_483
S 863 6 1 0 0 6 1 854 6633 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_486
S 864 23 5 0 0 8 866 624 6641 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 fftfreq
S 865 6 3 1 0 6 1 864 6649 800004 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 n
S 866 14 5 0 0 189 1 864 6641 204 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 12 1 0 0 867 0 0 0 0 0 0 0 0 0 218 0 624 0 0 0 0 fftfreq
F 866 1 865
S 867 7 3 0 0 189 1 864 6641 800204 1003000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 fftfreq
S 868 6 1 0 0 6 1 864 6651 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_252
S 869 23 5 0 0 8 871 624 6659 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 quad_shift
S 870 7 3 0 0 192 1 869 6670 20400004 10003000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 array
S 871 14 5 0 0 195 1 869 6659 20000204 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 14 1 0 0 872 0 0 0 0 0 0 0 0 0 235 0 624 0 0 0 0 quad_shift
F 871 1 870
S 872 7 3 0 0 195 1 869 6659 800204 1003000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 quad_shift
S 873 6 1 0 0 6 1 869 6498 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_0_1
S 874 6 1 0 0 6 1 869 6506 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_2_1
S 875 6 1 0 0 6 1 869 6514 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_3_1
S 876 6 1 0 0 6 1 869 6522 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_5_1
S 877 6 1 0 0 6 1 869 6530 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_b_6_1
S 878 6 1 0 0 6 1 869 6676 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_262
S 879 6 1 0 0 6 1 869 6684 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_265
S 880 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 881 6 1 0 0 6 1 869 6692 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_270
S 882 6 1 0 0 6 1 869 6700 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_274
S 883 6 1 0 0 6 1 869 6708 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_277
S 884 6 1 0 0 6 1 869 6716 40800006 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_279
A 231 1 0 0 0 6 841 0 0 0 0 0 0 0 0 0 0 0 0 0
A 232 1 0 0 0 6 837 0 0 0 0 0 0 0 0 0 0 0 0 0
A 233 1 0 0 0 6 842 0 0 0 0 0 0 0 0 0 0 0 0 0
A 234 1 0 0 0 6 839 0 0 0 0 0 0 0 0 0 0 0 0 0
A 235 1 0 0 0 6 838 0 0 0 0 0 0 0 0 0 0 0 0 0
A 236 1 0 0 0 6 843 0 0 0 0 0 0 0 0 0 0 0 0 0
A 237 1 0 0 0 6 840 0 0 0 0 0 0 0 0 0 0 0 0 0
A 238 1 0 0 0 6 851 0 0 0 0 0 0 0 0 0 0 0 0 0
A 239 1 0 0 9 6 847 0 0 0 0 0 0 0 0 0 0 0 0 0
A 240 1 0 0 0 6 852 0 0 0 0 0 0 0 0 0 0 0 0 0
A 241 1 0 0 0 6 849 0 0 0 0 0 0 0 0 0 0 0 0 0
A 242 1 0 0 0 6 848 0 0 0 0 0 0 0 0 0 0 0 0 0
A 243 1 0 0 0 6 853 0 0 0 0 0 0 0 0 0 0 0 0 0
A 244 1 0 0 0 6 850 0 0 0 0 0 0 0 0 0 0 0 0 0
A 245 1 0 0 0 6 861 0 0 0 0 0 0 0 0 0 0 0 0 0
A 246 1 0 0 0 6 857 0 0 0 0 0 0 0 0 0 0 0 0 0
A 247 1 0 0 97 6 862 0 0 0 0 0 0 0 0 0 0 0 0 0
A 248 1 0 0 0 6 859 0 0 0 0 0 0 0 0 0 0 0 0 0
A 249 1 0 0 0 6 858 0 0 0 0 0 0 0 0 0 0 0 0 0
A 250 1 0 0 0 6 863 0 0 0 0 0 0 0 0 0 0 0 0 0
A 251 1 0 0 0 6 860 0 0 0 0 0 0 0 0 0 0 0 0 0
A 252 1 0 0 0 6 865 0 0 0 0 0 0 0 0 0 0 0 0 0
A 253 1 0 0 0 6 868 0 0 0 0 0 0 0 0 0 0 0 0 0
A 254 1 0 0 0 6 877 0 0 0 0 0 0 0 0 0 0 0 0 0
A 255 1 0 0 0 6 873 0 0 0 0 0 0 0 0 0 0 0 0 0
A 256 1 0 0 0 6 878 0 0 0 0 0 0 0 0 0 0 0 0 0
A 257 1 0 0 0 6 875 0 0 0 0 0 0 0 0 0 0 0 0 0
A 258 1 0 0 0 6 874 0 0 0 0 0 0 0 0 0 0 0 0 0
A 259 1 0 0 0 6 879 0 0 0 0 0 0 0 0 0 0 0 0 0
A 260 1 0 0 0 6 876 0 0 0 0 0 0 0 0 0 0 0 0 0
A 261 1 0 0 0 6 884 0 0 0 0 0 0 0 0 0 0 0 0 0
A 262 1 0 0 0 0 426 0 0 0 0 0 0 0 0 0 0 0 0 0
A 263 1 0 7 0 192 870 0 0 0 0 0 0 0 0 0 0 0 0 0
A 264 14 0 0 0 6 262 0 0 0 0 0 0 243 2 1 0 0 0 0
W 2 263 3
A 265 1 0 0 0 6 881 0 0 0 0 0 0 0 0 0 0 0 0 0
A 266 2 0 0 0 6 880 0 0 0 266 0 0 0 0 0 0 0 0 0
A 267 14 0 0 0 6 262 0 0 0 0 0 0 243 2 4 0 0 0 0
W 2 263 266
A 268 1 0 0 0 6 882 0 0 0 0 0 0 0 0 0 0 0 0 0
A 269 1 0 0 0 6 883 0 0 0 0 0 0 0 0 0 0 0 0 0
Z
Z
