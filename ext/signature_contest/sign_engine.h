/// @file
/// @brief Определение интерфейса модулей, подаваемых на конкурс.
/// @copyright InfoTeCS. All rights reserved.

#ifndef CONTEST_SIGN_ENGINE_H
#define CONTEST_SIGN_ENGINE_H

#ifdef __cplusplus
extern "C" {
#endif //__cplusplus

typedef enum
{
     kStatusOk,
     kStatusWrongSignature,
     kStatusBadInput,
     kStatusInternalError
} Gost12S512Status;

/// @brief Инициализация модуля.
/// @return kStatusOk В случае успешного завершения.
/// @return kStatusInternalError В случае ошибки
Gost12S512Status Gost12S512Init( );

/// @brief Вычисление подписи сообщения для алгоритма хеширования Стрибог-512.
/// Все массивы, представляющие длинные целые числа, подаются на вход в формате LE. Все указатели выровнены.
/// @param[in] privateKey Закрытый ключ подписи.
/// @param[in] rand Случайное число.
/// @param[in] hash Подписываемый хеш.
/// @param[out] signature Сгенерированная подпись сообщения.
/// @return kStatusOk В случае успешного завершения.
/// @return kStatusBadInput Некорректные входные данные.
/// @return kStatusInternalError В остальных случаях
Gost12S512Status Gost12S512Sign( const char* privateKey,
                                 const char* rand,
                                 const char* hash,
                                 char* signature );



/// @brief Проверка корректности подписи сообщения для алгоритма хеширования Стрибог-512.
/// Все массивы, представляющие длинные целые числа, подаются на вход в формате LE. Все указатели выровнены.
/// @param[in] publicKeyX Точка X открытого ключа подписи.
/// @param[in] publicKeyУ Точка Y открытого ключа подписи.
/// @param[in] hash Хеш подписанного сообщения.
/// @param[in] signature Подпись.
/// @return kStatusOk В случае успешного завершения.
/// @return kStatusWrongSignature Подпись не соответствует исходному сообщению.
/// @return kStatusBadInput Некорректные входные данные.
/// @return kStatusInternalError В остальных случаях.
Gost12S512Status Gost12S512Verify( const char* publicKeyX,
                                   const char* publicKeyY,
                                   const char* hash,
                                   const char* signature );


#ifdef __cplusplus
}
#endif //__cplusplus

#endif //CONTEST_SIGN_ENGINE_H
