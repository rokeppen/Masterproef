function sol=euler_sol(t)
[SN,CN,DN] = ellipj(t,0.51); 
sol=[sqrt(1.51)*SN,CN,DN];
end