N = 11

print(f"Dimensions of M_k(Gamma_0({N})):\n")
print(f"{'Weight k':<10} {'dim M_k':<10}")
print("-" * 20)

for k in range(2, 21):
    dim = ModularForms(Gamma0(N), k).dimension()
    print(f"{k:<10} {dim:<10}")

print("\n" + "="*40 + "\n")

print("Detailed examples:\n")

for k in [2, 4, 6, 12]:
    dim_cusp = CuspForms(Gamma0(N), k).dimension()
    dim_modular = ModularForms(Gamma0(N), k).dimension()
    dim_eis = EisensteinForms(Gamma0(N), k).dimension()
    
    print(f"Weight k = {k}:")
    print(f"  dim M_k(Gamma_0({N})) = {dim_modular}")
    print(f"  dim S_k(Gamma_0({N})) = {dim_cusp}")
    print(f"  dim E_k(Gamma_0({N})) = {dim_eis}")
    print()

print("="*40)
print("\nFor weight 2, the space S_2(Gamma_0(11)):\n")

S = CuspForms(Gamma0(11), 2)
print(f"Dimension: {S.dimension()}")
print(f"q-expansion basis (first few terms):")
for i, f in enumerate(S.basis()):
    print(f"  f_{i+1} = {f.q_expansion(10)}")
