Input notes: 
Ran at L=1000 (up from L=500) whilst keeping kappa=1e-4. Also drastucally reduced the number of rho points on the base search.

Result notes:
Still no more steps found for rho>2.5, this must mean a smaller kappa is also needed at this higher L=1000, so will try kappa=1e-5 next at same params.

Reduced rho points works great, still found exactly the same steps up untul rho=2.5. Though, query, could the algorithm be missing steps above due to poor rho reolution? I imagine not, as the idea is to get high mu resolution at an arbitrarily higher rho, then check if there is a number of boundaries difference and check intermediate rho's until it satisfied (up to some depth limits).

Note, the lower rho resolutin does seem to have decreased the abaility of the algotihm to track boundaried correclty, there is much more crossing/jumping between visually clear boundaries in this run. Seems to leave many `ambiguous' and struggle to connect upper boundary of trivial regions more than lower.