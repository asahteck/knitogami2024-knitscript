//
//  ArgumentParser.hpp
//  Elasticity
//
//  Created by Wim van Rees on 8/25/16.
//  Copyright © 2016 Wim van Rees. All rights reserved.
//

#ifndef ArgumentParser_hpp
#define ArgumentParser_hpp


#include <iostream>
#include <cstring>
#include <map>
#include <string>
#include <sstream>
#include <ctype.h>

// ==== helper methods : convert string to type (and specializations for string and bool) ===== //
template<typename T>
inline T convertStrToType( const std::string & s )
{
    std::stringstream convert(s);

    T value;
    if(convert >> value)
        return value;
    else
    {
        std::cout << "argument value for type is invalid : \t " << s << std::endl;
        std::exit(1);
    }
    return value;
}

template<>
inline bool convertStrToType( const std::string & s )
{
    std::stringstream convert(s);

    bool value;
    if(convert >> std::boolalpha >> value)
        return value;
    else
    {
        std::cout << "argument value for bool type is invalid : \t " << s << std::endl;
        std::cout << "use true/false only" << std::endl;
        std::exit(1);
    }
    return value;
}

template<>
inline std::string convertStrToType(const std::string & s)
{
    return s;
}

// ==== helper methods : convert type to string (and specializations for string and bool) ===== //
template<typename T>
inline std::string convertTypeToStr(const T & t)
{
    return std::to_string(t);
}

template<>
inline std::string convertTypeToStr(const bool & t)
{
    return (t ? "true" : "false");
}

template<>
inline std::string convertTypeToStr(const std::string & t)
{
    return t;
}

// ==== base class Parser ==== //
class Parser
{
protected:
    std::map<std::string,std::string> arguments_map;
    bool verbose, saveDefaults;

    template<typename T>
    T parse_(const std::string argkey, const bool strict, const T defval = T()) //
    {
        std::map<std::string,std::string>::const_iterator it = arguments_map.find(argkey);

        if (it == arguments_map.end())
        {
            if (verbose) printf("%s is empty\n", argkey.data());
            if(strict)
            {
                helpers::catastrophe("Mandatory command line option is not given. Argument name: "+argkey, __FILE__, __LINE__);
            }
            if(saveDefaults)
            {
                FILE * f = getFileHandle();
                if(f!=nullptr)
                {
                    fprintf(f, "%s %s ", argkey.data(), convertTypeToStr(defval).data());
                    fclose(f);
                }

                // add to arguments map
                arguments_map[argkey] = convertTypeToStr(defval).data();
            }

            return defval;
        }
        else
        {
            if (verbose) printf("Found the value for key %s as %s\n", argkey.data(), it->second.data());
            return convertStrToType<T>(it->second);
        }
    }

    FILE * getFileHandle() const
    {
        const std::string filepath = "argumentparser.log";
        FILE * f = fopen(filepath.c_str(), "a");
        if (f == nullptr)
        {
            printf("Can not open file %s.\n", filepath.data());
            return nullptr;
        }
        return f;
    }

public:

    template<typename T>
    T parse(const std::string arg) //
    {
        return parse_<T>(arg, true);
    }

    template<typename T>
    T parse(const std::string arg, const T defval) //
    {
        return parse_<T>(arg, false, defval);
    }

    Parser(): arguments_map(), verbose(false), saveDefaults(false)
    {
    }

    void set_verbosity(const bool verbosity)
    {
        verbose = verbosity;
    }

    void save_defaults()
    {
        saveDefaults = true;
    }

    void print_options() const
    {
        for(const auto & it : arguments_map)
            printf("%s %s ",it.first.data(), it.second.data());
        printf("\n");
    }

    std::string options_string() const
    {
        std::stringstream ss;
        for(const auto & it : arguments_map) {
            ss << it.first.data() << " " << it.second.data() << " ";
        }
        return ss.str();
    }

    void save_options()
    {
        FILE * f = getFileHandle();
        if(f==nullptr) return;
        for(const auto & it : arguments_map) {
            fprintf(f, "%s %s ", it.first.data(), it.second.data());
        }
        fclose(f);
    }

    void initialize() const
    {
        FILE * f = getFileHandle();
        if(f==nullptr) return;
        fprintf(f, "\n\n");
        fclose(f);
    }

    void finalize() const
    {
        FILE * f = getFileHandle();
        if(f==nullptr) return;
        fprintf(f, "\n");
        fclose(f);
    }
};


class ArgumentParser
{
protected:
    std::vector< std::pair<std::string, std::string> > arguments_map;

    bool saveDefaults;

    const int parse_argc;
    const char** parse_argv;

    template<typename T>
    T parse_(const std::string argkey, const bool strict, const T defval = T())
    {
        // std::map<std::string,std::string>::const_iterator it = arguments_map.find(argkey);
        auto result = std::find_if(arguments_map.begin(), arguments_map.end(), [argkey](std::pair<std::string, std::string> pair){ return pair.first == argkey; });

        if (result == arguments_map.end()) {
            if (strict) {
                helpers::catastrophe("Mandatory command line option is not given. Argument name: "+argkey, __FILE__, __LINE__);
            }

            if (saveDefaults) {
                FILE * f = getFileHandle();
                if(f!=nullptr)
                {
                    fprintf(f, "%s %s ", argkey.data(), convertTypeToStr(defval).data());
                    fclose(f);
                }

                // add to arguments map
                arguments_map.push_back(std::make_pair(argkey, convertTypeToStr(defval).data()));
            }

            return defval;
        }
        else {
            return convertStrToType<T>(result->second);
        }
    }

    template<typename T>
    std::vector<T> parseAll_(const std::string argkey)
    {
        auto result = std::find_if(arguments_map.begin(), arguments_map.end(), [argkey](std::pair<std::string, std::string> pair){ return pair.first == argkey; });
        std::vector<T> vals;

        while (result != arguments_map.end()) {
            vals.push_back(convertStrToType<T>(result->second));
            result = std::find_if(++result, arguments_map.end(), [argkey](std::pair<std::string, std::string> pair){ return pair.first == argkey; });
        }

        return vals;
    }

    FILE * getFileHandle() const
    {
        const std::string filepath = "argumentparser.log";
        FILE * f = fopen(filepath.c_str(), "a");
        if (f == nullptr)
        {
            printf("Can not open file %s.\n", filepath.data());
            return nullptr;
        }
        return f;
    }

public:

    ArgumentParser(const int argc, const char ** argv) : parse_argc(argc), parse_argv(argv)
    {
        // loop over the arguments to build up arguments_map
        for (int i=1; i<argc; i++) {
            if (i == 1 and argv[i][0] != '-') {
                helpers::catastrophe("found an unexpected command-line entry : "+std::string(argv[i]), __FILE__, __LINE__);
            }

            if (i == argc-1) { // last field given with no value
                arguments_map.push_back(std::make_pair(argv[i], "true"));
                break;
            }

            if ( (argv[i+1][0] == '-') and not std::isdigit(argv[i+1][1]) ) { // no value given
                arguments_map.push_back(std::make_pair(argv[i], "true"));
                continue;
            }

            int j = i + 1;
            std::string values = "";
            while (j < argc) { // get value(s), at least one, for this argument
                values.append((j == i + 1 ? std::string("") : std::string(" ")) + std::string(argv[j]));

                if ( j == argc - 1 or ((argv[j+1][0] == '-') and not std::isdigit(argv[j+1][1])) ) { // next value is a new argument
                    break;
                }

                j++;
            }

            arguments_map.push_back(std::make_pair(argv[i], values));

            i = j;
        }

        std::cout << "arguments:" << std::endl;
        for (auto x : arguments_map) {
            std::cout << "\t" << x.first << ", " << x.second << std::endl;
        }
    }

    int getargc() const
    {
        return parse_argc;
    }

    const char** getargv() const
    {
        return parse_argv;
    }

    template<typename T>
    T add(const std::string arg, const T val)
    {
        arguments_map.push_back(std::make_pair(arg, convertTypeToStr(val).data()));
        return val;
    }

    template<typename T>
    T parse(const std::string arg)
    {
        return parse_<T>(arg, true);
    }

    template<typename T>
    T parse(const std::string arg, const T defval)
    {
        return parse_<T>(arg, false, defval);
    }

    template<typename T>
    std::vector<T> parseAll(const std::string arg)
    {
        return parseAll_<T>(arg);
    }

    void save_defaults()
    {
        saveDefaults = true;
    }

    void print_options() const
    {
        for(const auto & it : arguments_map)
            printf("%s %s ",it.first.data(), it.second.data());
        printf("\n");
    }

    std::string options_string() const
    {
        std::stringstream ss;
        for(const auto & it : arguments_map) {
            ss << it.first.data() << " " << it.second.data() << " ";
        }
        return ss.str();
    }

    void save_options()
    {
        FILE * f = getFileHandle();
        if(f==nullptr) return;
        for(const auto & it : arguments_map) {
            fprintf(f, "%s %s ", it.first.data(), it.second.data());
        }
        fclose(f);
    }

    void initialize() const
    {
        FILE * f = getFileHandle();
        if(f==nullptr) return;
        fprintf(f, "\n");
        fclose(f);
    }

    void finalize() const
    {
        // FILE * f = getFileHandle();
        // if(f==nullptr) return;
        // fprintf(f, "\n");
        // fclose(f);
    }
};

#endif /* ArgumentParser_hpp */
