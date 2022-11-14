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
Integer Variable
    u(i,k)  'Additional variable for subtour elimination';
    
Variables
    Z      'target value as total sum of transport costs';

Equations
    target_function                     'target function'
    demand_capa_constr(k)               'demands assigned to a transport vehicle correspond to max. vehicle capacity'
    depot_constr(i)                     'tours include depot as start and end'
    cust_vehicle_constr(i)              'customer is assigned to exactly one transport vehicle'
    node_after_constr(i,k)              'each node has one node afterwards'
    node_before_const(j,k)              'each node has one node before'
    subtour_elimination_constr(i,j,k)   'avoid subtours without depot'
    lower_bound_subtour_constr(i,k)     'lower bound to avoid subtours without depot'   
    upper_bound_subtour_constr(i,k)     'upper bound to avoid subtours without depot';
      
target_function..                                                                                       Z =e= sum((i,j,k), c(i,j)*x(i,j,k));
demand_capa_constr(k)..                                                                                 sum(i, b(i)*y(i,k)) =l= a(k);
depot_constr(i)$((ord(i)=1))..                                                                          sum(k, y(i,k)) =e= card(k);
cust_vehicle_constr(i)$((ord(i)<>1))..                                                                  sum(k, y(i,k)) =e= 1;
node_after_constr(i,k)..                                                                                sum(j, x(i,j,k)) =e= y(i,k);
node_before_const(j,k)..                                                                                sum(i, x(i,j,k)) =e= y(j,k);
subtour_elimination_constr(i,j,k)$(not sameAs(i,j)and(ord(i)<>1)and(ord(j)<>1)and(b(i)+b(j)<=a(k)))..   u(i,k)-u(j,k)+a(k)*x(i,j,k) =l= a(k)-b(j);             
lower_bound_subtour_constr(i,k)..                                                                       b(i)=l=u(i,k);
upper_bound_subtour_constr(i,k)..                                                                       u(i,k)=l=a(k) ;

Model CVRP /all/;

Solve CVRP using MIP minimizing Z;