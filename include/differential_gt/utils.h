
namespace std
{
inline std::string to_string ( const std::vector<std::string>& vals )
{
    std::string ret = "< ";
    for ( auto const & val : vals ) ret += val + ", ";
    ret += " >";
    return ret;
}

inline std::string to_string ( const std::string& val )
{
    return val;
}

inline std::string to_string ( const bool& val )
{
    return val ? "TRUE" : "FALSE";
}

}


#define GET_AND_RETURN( nh, param, value )\
  if (!nh.getParam(param,value) )\
  {\
    ROS_ERROR("The param '%s/%s' is not defined", nh.getNamespace().c_str(), std::string( param ).c_str() );\
    return false;\
  }



#define GET_AND_DEFAULT( nh, param, value, def )\
  if (!nh.getParam(param,value) )\
  {\
    ROS_WARN("The param '%s/%s' is not defined", nh.getNamespace().c_str(), std::string( param ).c_str() );\
    ROS_WARN("Default value '%s' superimposed. ", std::to_string( def ).c_str() );\
    value=def;\
  }

#define GET_PARAM_VECTOR_AND_RETURN(nh, P, X , N, DEQ )\
  if (!nh.getParam( std::string(P).c_str(), X))\
  {\
    ROS_FATAL_STREAM("[ " << nh.getNamespace().c_str() << "] Parameter '"<<  P <<"' does not exist");\
    ROS_FATAL_STREAM("[ " << nh.getNamespace().c_str() << "] ERROR DURING INITIALIZATION. ABORT.");\
    return false;\
  }\
  if( X.size() != N )\
  {\
    ROS_FATAL_STREAM("[ " << nh.getNamespace().c_str() << "] The size '"<< X.size() <<"' of the param '" << P << "' does not match with the foreseen dimension '"<< N <<"'");\
    ROS_FATAL_STREAM("[ " << nh.getNamespace().c_str() << "] ERROR DURING INITIALIZATION. ABORT.");\
    return false;\
  }\
  for( size_t i=0; i<X.size(); i++)\
  {\
    if( std::string(DEQ)==std::string("<"))\
    {\
      if( X.at(i) < 0) \
      {\
        ROS_FATAL_STREAM("[ " << nh.getNamespace().c_str() << "] Parameter '"<<  P <<"' has negative values. Abort");\
        return false;\
      }\
    }\
    else if( std::string(DEQ)==std::string("<="))\
    {\
      if( X.at(i) <= 0) \
      {\
        ROS_FATAL_STREAM("[ " << nh.getNamespace().c_str() << "] Parameter '"<<  P <<"' has negative or null values. Abort");\
        return false;\
      }\
    }\
  }\

namespace eu = eigen_utils;
namespace ect = eigen_control_toolbox;


static const char* DEFAULT      = "\033[0m";
static const char* RESET        = "\033[0m";
static const char* BLACK        = "\033[30m";
static const char* RED          = "\033[31m";
static const char* GREEN        = "\033[32m";
static const char* YELLOW       = "\033[33m";
static const char* BLUE         = "\033[34m";
static const char* MAGENTA      = "\033[35m";
static const char* CYAN         = "\033[36m";
static const char* WHITE        = "\033[37m";
static const char* BOLDBLACK    = "\033[1m\033[30m";
static const char* BOLDRED      = "\033[1m\033[31m";
static const char* BOLDGREEN    = "\033[1m\033[32m";
static const char* BOLDYELLOW   = "\033[1m\033[33m";
static const char* BOLDBLUE     = "\033[1m\033[34m";
static const char* BOLDMAGENTA  = "\033[1m\033[35m";
static const char* BOLDCYAN     = "\033[1m\033[36m";
static const char* BOLDWHITE    = "\033[1m\033[37m";
