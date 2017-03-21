V29 :0x4 cuda_setup
18 mod_cuda_setup.f90 S624 0
03/07/2017  12:04:36
use iso_c_binding public 0 indirect
use pgi_acc_common public 0 indirect
use cudafor_lib public 0 indirect
use cudafor public 0 direct
use m_user_input private
enduse
D 56 24 646 8 645 7
D 62 24 648 8 647 7
D 74 24 646 8 645 7
D 92 24 721 8 720 7
S 624 24 0 0 0 8 1 0 5011 10005 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 cuda_setup
S 627 23 0 0 0 8 10676 624 5043 4 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 get_input
R 645 25 6 iso_c_binding c_ptr
R 646 5 7 iso_c_binding val c_ptr
R 647 25 8 iso_c_binding c_funptr
R 648 5 9 iso_c_binding val c_funptr
R 681 6 42 iso_c_binding c_null_ptr$ac
R 683 6 44 iso_c_binding c_null_funptr$ac
R 684 26 45 iso_c_binding ==
R 686 26 47 iso_c_binding !=
R 720 25 5 pgi_acc_common c_devptr
R 721 5 6 pgi_acc_common cptr c_devptr
R 723 6 8 pgi_acc_common c_null_devptr$ac
R 727 26 12 pgi_acc_common =
R 10676 19 1 m_user_input get_input
S 10777 6 4 0 0 7 1 624 69948 4 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 10778 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 usable_gpu_memory
S 10778 11 0 0 0 8 10697 624 69966 40800000 805000 A 0 0 0 0 B 0 0 0 0 0 8 0 0 10777 10777 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _cuda_setup$2
S 10779 23 5 0 0 0 10781 624 69980 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 display_device_properties
S 10780 1 3 0 0 6 1 10779 6044 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 i
S 10781 14 5 0 0 0 1 10779 69980 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 4790 1 0 0 0 0 0 0 0 0 0 0 0 0 12 0 624 0 0 0 0 display_device_properties
F 10781 1 10780
S 10782 23 5 0 0 0 10783 624 70006 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 setup_gpu
S 10783 14 5 0 0 0 1 10782 70006 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 4792 0 0 0 0 0 0 0 0 0 0 0 0 0 44 0 624 0 0 0 0 setup_gpu
F 10783 0
S 10784 23 5 0 0 0 10787 624 70016 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 gpu_memory_message
S 10785 1 3 1 0 9 1 10784 70035 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 required_memory
S 10786 1 3 2 0 16 1 10784 70051 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 on_the_fly
S 10787 14 5 0 0 0 1 10784 70016 0 400000 A 0 0 0 0 B 0 0 0 0 0 0 0 4793 2 0 0 0 0 0 0 0 0 0 0 0 0 126 0 624 0 0 0 0 gpu_memory_message
F 10787 2 10785 10786
A 67 1 0 0 0 56 681 0 0 0 0 0 0 0 0 0 0 0 0 0
A 70 1 0 0 0 62 683 0 0 0 0 0 0 0 0 0 0 0 0 0
A 86 1 0 0 0 92 723 0 0 0 0 0 0 0 0 0 0 0 0 0
Z
J 149 1 1
V 67 56 7 0
S 0 56 0 0 0
A 0 6 0 0 1 2 0
J 150 1 1
V 70 62 7 0
S 0 62 0 0 0
A 0 6 0 0 1 2 0
J 31 1 1
V 86 92 7 0
S 0 92 0 0 0
A 0 74 0 0 1 67 0
Z
