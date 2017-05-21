        .intel_syntax noprefix
        .section .text
        .global rsqrt
        .p2align 4,,15
rsqrt:
        mov     ecx, 64
        vzeroall
rsqrt_loop:
        vmovapd ymm12, [rdi]
#        vcvtpd2ps xmm1, ymm0
#        vrsqrtps xmm2, xmm1
#        vcvtps2pd ymm3, xmm2
#        vbroadcastsd ymm4, .onehalf[rip]
#        vmulpd  ymm5, ymm0, .neghalf[rip]
#
#        vmulpd  ymm6, ymm5, ymm3
#        vfmadd213pd ymm6, ymm3, ymm4
#        vmulpd  ymm7, ymm6, ymm3
#        vmulpd  ymm8, ymm5, ymm7
#        vfmadd213pd ymm8, ymm7, ymm4
#        vmulpd  ymm0, ymm8, ymm7

        vcvtpd2ps xmm13, ymm12
        vrsqrtps xmm13, xmm13
        vcvtps2pd ymm13, xmm13
        vbroadcastsd ymm7, .onehalf[rip]
        vmulpd  ymm12, ymm12, .neghalf[rip]

        vmulpd  ymm14, ymm13, ymm12
        vfmadd213pd ymm14, ymm13, ymm7
        vmulpd  ymm13, ymm14, ymm13
        vmulpd  ymm14, ymm13, ymm12
        vfmadd213pd ymm14, ymm13, ymm7
        vmulpd  ymm13, ymm14, ymm13

        vmovapd [rdi], ymm13
        add     rdi, 32
        loop    rsqrt_loop
        ret
        .section .rodata
        .align  8
.onehalf:
        .double 1.5
        .align  32
.neghalf:
        .double -0.5
        .double -0.5
        .double -0.5
        .double -0.5
