$Set NODES 51
Set
    i   "nodos"
        /depot,N1*N%NODES%/
    s(i)  "depot node"
        /depot/
    t(i)  "nodos fundaciones"
        /N1*N%NODES%/
    k   "CARGA"
        /T1*T8/
    ;
alias(j,i);

Set r(j) 'customer nodes'  /N1*N%NODES%/

Table c(i,j)
$ondelim
$include distancias.csv
$offdelim
;
Display c;


Parameters b(i) 'Demanda de mercados por fundación'
/
$ondelim
$include familias_atendidas.csv
$offdelim
/
;
Display b;
Parameters a(k) 'Capacidad de transporte del vehiculo'
/
$ondelim
$include camiones.csv
$offdelim
/
;
Display a;
Scalar M 'cantidad de vehiculos de transporte' /8/;
Binary variable
x(i,j,k) '1 Si el camion k viajo desde el nodo i al nodo j, sino 0'
y(i,k)  '1 SI el cliente i fue completamente avastesido por el camion k, sino o'

;
Variables
u(t)
v(r)
z;


Equations
   target_function                            'target function'
   demand_capa_constraint(k)                  'demands assigned to a transport vehicle correspond to max. vehicle capacity'
   depot_constraint(s)                        'tours include depot as start and end'
   customer_vehicle_constraint(t)             'customer is assigned to exactly one transport vehicle'
   node_after_constraint(i,k)                 'each node has one node afterwards'
   node_before_constraint(j,k)                'each node has one node before'
   subtour_elimination_constraint(r,t,i,j,k)  'avoid subtours without depot'
   subtour_elimination_constraint2(i,t)       'constraint for subtour elimination'
   subtour_elimination_constraint3(t,k)       'constraint for subtour elimination';
   
target_function..                           Z =e= sum((i,j,k), c(i,j)*x(i,j,k));
demand_capa_constraint(k)..                 sum(i, b(i)*y(i,k)) =l= a(k);
depot_constraint(s)..                       M =e= sum(k, y(s,k));
customer_vehicle_constraint(t)..            sum(k, y(t,k)) =e= 1;
node_after_constraint(i,k)..                sum(j, x(i,j,k)) =e= y(i,k);
node_before_constraint(j,k)..               sum(i, x(i,j,k)) =e= y(j,k);
subtour_elimination_constraint(r,t,i,j,k).. v(r) - u(t) =g= b(j) - a(k) * (1-x(i,j,k));
subtour_elimination_constraint2(i,t)..      b(i)=l= u(t);
subtour_elimination_constraint3(t,k)..      u(t)=l= a(k);

Model CVRP /all/;

Solve CVRP using MIP minimizing Z
