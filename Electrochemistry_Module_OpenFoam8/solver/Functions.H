/*namespace electro {

template <class T>
inline T call_doubleWellPot(T& eta, T& phi) {
	const auto R = 3.13;
	return 2*phi*eta * (1-eta)*(1-2*eta);
}



};

electro::call_doubleWellPot();*/
template <class T>
inline T call_doubleWellPot(T& eta) {
	//const auto W = 3.13;
	return 2*eta*(1-eta)*(1-2*eta);
}
template <class T>
inline T call_InterPotFn(T& eta) {
	//const auto W = 3.13;
	return 30*eta*eta*(1-eta)*(1-eta);
}
template <class T>
inline T call_Conc_Li(T& mu, T& eta) {
	//const auto W = 3.13;
	return (1-(eta*eta*eta*(6*eta*eta - 15*eta +10)))*exp(mu);
}
