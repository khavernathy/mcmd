/****************************************************************************
** Meta object code from reading C++ file 'fileio.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.9.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../gui/fileio.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'fileio.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.9.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_FileIO_t {
    QByteArrayData data[21];
    char stringdata0[189];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_FileIO_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_FileIO_t qt_meta_stringdata_FileIO = {
    {
QT_MOC_LITERAL(0, 0, 6), // "FileIO"
QT_MOC_LITERAL(1, 7, 13), // "sourceChanged"
QT_MOC_LITERAL(2, 21, 0), // ""
QT_MOC_LITERAL(3, 22, 6), // "source"
QT_MOC_LITERAL(4, 29, 15), // "grSourceChanged"
QT_MOC_LITERAL(5, 45, 8), // "grSource"
QT_MOC_LITERAL(6, 54, 16), // "linecountChanged"
QT_MOC_LITERAL(7, 71, 4), // "int&"
QT_MOC_LITERAL(8, 76, 9), // "linecount"
QT_MOC_LITERAL(9, 86, 14), // "homeDirChanged"
QT_MOC_LITERAL(10, 101, 7), // "homeDir"
QT_MOC_LITERAL(11, 109, 5), // "error"
QT_MOC_LITERAL(12, 115, 3), // "msg"
QT_MOC_LITERAL(13, 119, 9), // "setSource"
QT_MOC_LITERAL(14, 129, 12), // "setLinecount"
QT_MOC_LITERAL(15, 142, 10), // "setHomeDir"
QT_MOC_LITERAL(16, 153, 11), // "setgrSource"
QT_MOC_LITERAL(17, 165, 4), // "read"
QT_MOC_LITERAL(18, 170, 7), // "read_gr"
QT_MOC_LITERAL(19, 178, 5), // "write"
QT_MOC_LITERAL(20, 184, 4) // "data"

    },
    "FileIO\0sourceChanged\0\0source\0"
    "grSourceChanged\0grSource\0linecountChanged\0"
    "int&\0linecount\0homeDirChanged\0homeDir\0"
    "error\0msg\0setSource\0setLinecount\0"
    "setHomeDir\0setgrSource\0read\0read_gr\0"
    "write\0data"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_FileIO[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      12,   14, // methods
       4,  106, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       5,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,   74,    2, 0x06 /* Public */,
       4,    1,   77,    2, 0x06 /* Public */,
       6,    1,   80,    2, 0x06 /* Public */,
       9,    1,   83,    2, 0x06 /* Public */,
      11,    1,   86,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
      13,    1,   89,    2, 0x0a /* Public */,
      14,    1,   92,    2, 0x0a /* Public */,
      15,    1,   95,    2, 0x0a /* Public */,
      16,    1,   98,    2, 0x0a /* Public */,

 // methods: name, argc, parameters, tag, flags
      17,    0,  101,    2, 0x02 /* Public */,
      18,    0,  102,    2, 0x02 /* Public */,
      19,    1,  103,    2, 0x02 /* Public */,

 // signals: parameters
    QMetaType::Void, QMetaType::QString,    3,
    QMetaType::Void, QMetaType::QString,    5,
    QMetaType::Void, 0x80000000 | 7,    8,
    QMetaType::Void, QMetaType::QString,   10,
    QMetaType::Void, QMetaType::QString,   12,

 // slots: parameters
    QMetaType::Void, QMetaType::QString,    3,
    QMetaType::Void, 0x80000000 | 7,    8,
    QMetaType::Void, QMetaType::QString,   10,
    QMetaType::Void, QMetaType::QString,    5,

 // methods: parameters
    QMetaType::QString,
    QMetaType::QString,
    QMetaType::Bool, QMetaType::QString,   20,

 // properties: name, type, flags
       3, QMetaType::QString, 0x00495103,
       5, QMetaType::QString, 0x00495003,
       8, QMetaType::Int, 0x00495103,
      10, QMetaType::QString, 0x00495103,

 // properties: notify_signal_id
       0,
       1,
       2,
       3,

       0        // eod
};

void FileIO::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        FileIO *_t = static_cast<FileIO *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->sourceChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 1: _t->grSourceChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 2: _t->linecountChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 3: _t->homeDirChanged((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 4: _t->error((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 5: _t->setSource((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 6: _t->setLinecount((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 7: _t->setHomeDir((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 8: _t->setgrSource((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 9: { QString _r = _t->read();
            if (_a[0]) *reinterpret_cast< QString*>(_a[0]) = std::move(_r); }  break;
        case 10: { QString _r = _t->read_gr();
            if (_a[0]) *reinterpret_cast< QString*>(_a[0]) = std::move(_r); }  break;
        case 11: { bool _r = _t->write((*reinterpret_cast< const QString(*)>(_a[1])));
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = std::move(_r); }  break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (FileIO::*_t)(const QString & );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&FileIO::sourceChanged)) {
                *result = 0;
                return;
            }
        }
        {
            typedef void (FileIO::*_t)(const QString & );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&FileIO::grSourceChanged)) {
                *result = 1;
                return;
            }
        }
        {
            typedef void (FileIO::*_t)(int & );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&FileIO::linecountChanged)) {
                *result = 2;
                return;
            }
        }
        {
            typedef void (FileIO::*_t)(const QString & );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&FileIO::homeDirChanged)) {
                *result = 3;
                return;
            }
        }
        {
            typedef void (FileIO::*_t)(const QString & );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&FileIO::error)) {
                *result = 4;
                return;
            }
        }
    }
#ifndef QT_NO_PROPERTIES
    else if (_c == QMetaObject::ReadProperty) {
        FileIO *_t = static_cast<FileIO *>(_o);
        Q_UNUSED(_t)
        void *_v = _a[0];
        switch (_id) {
        case 0: *reinterpret_cast< QString*>(_v) = _t->source(); break;
        case 1: *reinterpret_cast< QString*>(_v) = _t->grSource(); break;
        case 2: *reinterpret_cast< int*>(_v) = _t->linecount(); break;
        case 3: *reinterpret_cast< QString*>(_v) = _t->homeDir(); break;
        default: break;
        }
    } else if (_c == QMetaObject::WriteProperty) {
        FileIO *_t = static_cast<FileIO *>(_o);
        Q_UNUSED(_t)
        void *_v = _a[0];
        switch (_id) {
        case 0: _t->setSource(*reinterpret_cast< QString*>(_v)); break;
        case 1: _t->setgrSource(*reinterpret_cast< QString*>(_v)); break;
        case 2: _t->setLinecount(*reinterpret_cast< int*>(_v)); break;
        case 3: _t->setHomeDir(*reinterpret_cast< QString*>(_v)); break;
        default: break;
        }
    } else if (_c == QMetaObject::ResetProperty) {
    }
#endif // QT_NO_PROPERTIES
}

const QMetaObject FileIO::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_FileIO.data,
      qt_meta_data_FileIO,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *FileIO::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *FileIO::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_FileIO.stringdata0))
        return static_cast<void*>(const_cast< FileIO*>(this));
    return QObject::qt_metacast(_clname);
}

int FileIO::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 12)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 12;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 12)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 12;
    }
#ifndef QT_NO_PROPERTIES
   else if (_c == QMetaObject::ReadProperty || _c == QMetaObject::WriteProperty
            || _c == QMetaObject::ResetProperty || _c == QMetaObject::RegisterPropertyMetaType) {
        qt_static_metacall(this, _c, _id, _a);
        _id -= 4;
    } else if (_c == QMetaObject::QueryPropertyDesignable) {
        _id -= 4;
    } else if (_c == QMetaObject::QueryPropertyScriptable) {
        _id -= 4;
    } else if (_c == QMetaObject::QueryPropertyStored) {
        _id -= 4;
    } else if (_c == QMetaObject::QueryPropertyEditable) {
        _id -= 4;
    } else if (_c == QMetaObject::QueryPropertyUser) {
        _id -= 4;
    }
#endif // QT_NO_PROPERTIES
    return _id;
}

// SIGNAL 0
void FileIO::sourceChanged(const QString & _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void FileIO::grSourceChanged(const QString & _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void FileIO::linecountChanged(int & _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void FileIO::homeDirChanged(const QString & _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 3, _a);
}

// SIGNAL 4
void FileIO::error(const QString & _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 4, _a);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
