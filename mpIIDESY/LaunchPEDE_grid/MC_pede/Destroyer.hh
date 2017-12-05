#ifndef DESTROYER_HH
#define DESTROYER_HH

template <class DOOMED>

class Destroyer {
public:
  Destroyer(DOOMED* = 0);
  ~Destroyer();

  void SetDoomed(DOOMED*);
private:
// Prevent users from making copies of a
// Destroyer to avoid double deletion:
  Destroyer(const Destroyer<DOOMED>&);
  Destroyer<DOOMED>& operator=(const Destroyer<DOOMED>&);
  DOOMED* _doomed;
};


//
// Implement template class in header
//

template <class DOOMED>
Destroyer<DOOMED>::Destroyer (DOOMED* d) {
  _doomed = d;
}

template <class DOOMED>
Destroyer<DOOMED>::~Destroyer () {
  delete _doomed;
}

template <class DOOMED>
void Destroyer<DOOMED>::SetDoomed (DOOMED* d) {
  _doomed = d;
}

//template class Destroyer<Logger>;
//template class Destroyer<CmdParser>;
//template class Destroyer<DBase>;

#endif
