#include "random.hh"

class RNG_music : public RNG_plugin{
public:
  explicit RNG_music( config_file& cf )
  : RNG_plugin( cf )
  { }
    
  ~RNG_music() { }
    
  bool is_multiscale() const
  {   return true;   } 
};


namespace{
  RNG_plugin_creator_concrete< RNG_music > creator("MUSIC");


}
