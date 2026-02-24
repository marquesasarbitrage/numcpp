#pragma once 
#include <string>
#include <stdexcept>

namespace numcpp {

    namespace errors {

        enum class ErrorCode {

            UnknownError = 0,
            InvalidInput = 1, 
            NumericalError = 2, 
            UnknownTypeFromEnum = 6
        };

        inline std::string getCodeAsString(ErrorCode code) {
            switch (code) {
                case ErrorCode::UnknownError: return "UnknownError";
                case ErrorCode::InvalidInput: return "InvalidInput";
                case ErrorCode::NumericalError: return "NumericalError";
                case ErrorCode::UnknownTypeFromEnum: return "UnknownTypeFromEnum";
                default: return "UnknownError";
            }
        }

        inline std::string getCodeMessage(ErrorCode code) {
            switch (code) {
                case ErrorCode::UnknownError: return "Unknown error occured";
                case ErrorCode::InvalidInput: return "The input in the object or method is not valid";
                case ErrorCode::NumericalError: return "There was a numerical error in the computations";
                case ErrorCode::UnknownTypeFromEnum: return "The type as enum object is unknown from the used method";
                default: return "UnknownError";
            }
        }

        class Error : public std::runtime_error {

            public: 
                explicit Error(ErrorCode code): std::runtime_error(generateErrorMessage(code)), code_(code) {}
                explicit Error(ErrorCode code, std::string detail): std::runtime_error(generateErrorMessage(code, detail)), code_(code) {}
                ~Error() = default; 

                ErrorCode getCode() const noexcept {return code_;}
            private:
                ErrorCode code_;
                static std::string generateErrorMessage(ErrorCode code) { return getCodeAsString(code) + "::" + std::to_string(static_cast<int>(code)) + " : " + getCodeMessage(code); }
                static std::string generateErrorMessage(ErrorCode code, std::string detail) {return generateErrorMessage(code) + "[ Details : " + detail + " ]";}
        
        };


        inline void throwUnknownError() { throw Error(ErrorCode::UnknownError);}
        inline void throwUnknownError(std::string details) { throw Error(ErrorCode::UnknownError,details);}

        inline void throwInvalidInput() { throw Error(ErrorCode::InvalidInput);}
        inline void throwInvalidInputr(std::string details) { throw Error(ErrorCode::InvalidInput,details);}

        inline void throwNumericalError() { throw Error(ErrorCode::NumericalError);}
        inline void throwNumericalError(std::string details) { throw Error(ErrorCode::NumericalError,details);}

        inline void throwUnknownTypeFromEnum() { throw Error(ErrorCode::UnknownTypeFromEnum);}
        inline void throwUnknownTypeFromEnum(std::string details) { throw Error(ErrorCode::UnknownTypeFromEnum,details);}


    }

    
    
};