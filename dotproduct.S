        .intel_syntax noprefix
        .section .text
        .p2align 4,,15
        .global dotproduct
dotproduct:
        sub     rsp, 8
        shr     ecx, 2
dotproduct_loop:
        vmovapd ymm0, [rdi]
        vmovapd ymm1, 32[rdi]
        vmovapd ymm2, 64[rdi]
        vmovapd ymm3, 96[rdi]
        vmovapd ymm4, [rsi]
        vmovapd ymm5, 32[rsi]
        vmovapd ymm6, 64[rsi]
        vmovapd ymm7, 96[rsi]
        vmulpd  ymm0, ymm0, ymm4
        vmulpd  ymm1, ymm1, ymm5
        vmulpd  ymm2, ymm2, ymm6
        vmulpd  ymm3, ymm3, ymm7
        vhaddpd ymm0, ymm0, ymm1
        vhaddpd ymm1, ymm2, ymm3
        vperm2f128 ymm2, ymm0, ymm1, 0x21
        vblendpd ymm3, ymm0, ymm1, 0b1100
        vaddpd  ymm0, ymm2, ymm3
        vmovapd [rdx], ymm0
        add     rdi, 128
        add     rsi, 128
        add     rdx, 32
        loop    dotproduct_loop
        add     rsp, 8
        ret
