// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_XML_H
#define HEMELB_IO_XML_H

#include <cstdlib>
#include <cstring>
#include <functional>
#include <iostream>
#include <sstream>
#include <limits>
#include <memory>
#include <optional>
#include <filesystem>

#include "Exception.h"
#include "util/traits.h"

namespace tinyxml2 {
    // Forward declare the TinyXML types needed.
    class XMLDocument;
    class XMLElement;
    class XMLAttribute;
    // Worth noting that TinyXML-2 only works with C-style null terminated strings.
}

namespace hemelb::io::xml
{
    // Forward declare
    class Document;
    class NamedChildIterator;
    class UnnamedChildIterator;
    struct NamedIterationRange;
    struct UnnamedIterationRange;
    struct NamedPopIterationRange;
    struct AttributeNameIterationRange;
    class DyingElement;

    // Represent an element in an XML file, or no element (if equal to
    // the result of Missing());
    //
    // Exposes a basic monadic interface.
    class Element
    {
        tinyxml2::XMLElement* el = nullptr;

    public:

        static Element Missing();

        // The default Element() is equal to Element::Missing()
        Element() = default;

        Element(tinyxml2::XMLElement* el);
        ~Element();

        /**
         * Get the name of the element
         * @return
         */
        [[nodiscard]] std::string_view GetName() const;

        /**
         * Get the line number
         * @return
         */
        [[nodiscard]] int GetLine() const;

        /**
         * Gets the first child element, whatever its name.
         * @return element or missing if none
         */
        [[nodiscard]] Element GetChildOrNull() const;

        /**
         * Gets the first child element with the specified name
         *
         * @param $name
         *   the name of the child element to return
         *
         * @return
         *   returns the child element if it was found or
         *   Element::Missing() if not
         */
        Element GetChildOrNull(char const* name) const;

        /**
         * Gets the first child element whatever its name or throw
         * ChildError if it does not exist.
         *
         * @param $name
         *   the name of the child element to return
         *
         * @return
         *   returns the child element
         */
        [[nodiscard]] Element GetChildOrThrow() const;

        /**
         * Gets the first child element with the specified name or throw
         * ChildError if it does not exist.
         *
         * @param $name
         *   the name of the child element to return
         *
         * @return
         *   returns the child element
         */
        Element GetChildOrThrow(char const* name) const;

        /**
         * Return iterator over children with the specified name.
         * @param name
         * @return
         */
        [[nodiscard]] NamedIterationRange Children(char const* name) const;
        [[nodiscard]] NamedIterationRange Children(std::string name) const;

        [[nodiscard]] NamedPopIterationRange PopChildren(char const* name) const;
        [[nodiscard]] NamedPopIterationRange PopChildren(std::string name) const;

        /**
         * Return iterator over all children.
         * @param name
         * @return
         */
        [[nodiscard]] UnnamedIterationRange Children() const;

        /**
         * Return the next sibling element with the specified name, if any
         * If no suitable sibling element exists then return Element::Missing()
         *
         * @param $name
         *   the name of the sibling element to get
         *
         * @return
         *   Next sibling or Element::Missing()
         */
        [[nodiscard]] Element NextSiblingOrNull() const;

        /**
         * Return the next sibling element with the specified name, if any
         * If no suitable sibling element exists then return Element::Missing()
         *
         * @param $name
         *   the name of the sibling element to get
         *
         * @return
         *   Next sibling or Element::Missing()
         */
        [[nodiscard]] Element NextSiblingOrNull(char const* name) const;

        /**
         * Return the parent element unless this is the root element.
         * If so then return Element::Missing()
         *
         * @return
         *   Parent or Element::Missing()
         */
        [[nodiscard]] Element GetParentOrNull() const;
        /**
         * Return the parent element unless this is the root element.
         * If so then throws ParentError
         *
         * @return
         *   Parent
         */
        [[nodiscard]] Element GetParentOrThrow() const;

        // True if this element has at least one attribute.
        [[nodiscard]] bool HasAttributes() const;

        [[nodiscard]] AttributeNameIterationRange Attributes() const;

        /**
         * Get the value (as a string) contained in the specified attribute.
         * If it does not exist, return std::nullopt
         * @param $name
         *   The name of the attribute to get
         * @return
         *   An optional string containing the attribute value (or nullopt
         *   on failure)
         */
        [[nodiscard]] std::optional<std::string_view> GetAttributeMaybe(char const* name) const;

        /**
         * Get the value (as a string) contained in the specified attribute.
         * If it does not exist, throws AttributeError
         * @param $name
         *   The name of the attribute to get
         * @return
         *   A reference to a string containing the attribute value
         */
        [[nodiscard]] std::string_view GetAttributeOrThrow(char const* name) const;

        /**
         * Get the value contained in the specified attribute. This function
         * will attempt to convert the string to the type of templace
         * argument using std::istream::operator>> (i.e. the standard library's
         * formatted input operator) so this can work for arbitrary types,
         * as long as you supply an implementation for:
         *   std::istream& operator>>(std::istream&, T&);
         *
         * If the attribute is missing, you will get std::nullopt.
         *
         * If the attribute exists but does not convert correctly to type T,
         * then this will throw std::stringstream::failure.
         *
         * If the attribute exists, and converts to type T but does not use
         * the whole attribute string, then this will throw ParseError.
         *
         * @param $name
         *   Attribute to read and convert
         * @param $out
         *   Variable in which to store the converted attribute.
         * @return
         *   bool indicating if the attribute was present.
         */
        template<class T>
        [[nodiscard]] std::optional<T> GetAttributeMaybe(char const* name) const;

        template<class T>
        [[nodiscard]] T GetAttributeOrThrow(char const* name) const;

        /**
         * Get the value contained in the specified attribute. This function
         * will attempt to convert the string to the type of the second
         * argument using std::istream::operator>> (i.e. the standard library's
         * formatted input operator) so this can work for arbitrary types,
         * as long as you supply an implementation for:
         *   std::istream& operator>>(std::istream&, T&);
         *
         * This also returns the attribute as a reference to the string value
         * of the attribute.
         *
         * If the attribute does not exists, throw AttributeError.
         *
         * If the attribute exists but does not convert correctly to type T,
         * then this will throw std::stringstream::failure.
         *
         * If the attribute exists, and converts to type T but does not use
         * the whole attribute string, then this will throw DeserialisationError.
         *
         * @param $name
         *   Attribute to read and convert
         * @param $out
         *   Variable in which to store the converted attribute.
         * @return
         *   A reference to a string containing the attribute value
         */
        template<class T>
        void GetAttributeOrThrow(char const* name, T& out) const;

        void DelAttribute(char const* name);

        // These "pop" variants act as the "get" versions above, but also
        // delete the attribute.
        [[nodiscard]] std::optional<std::string> PopAttributeMaybe(char const* name);
        [[nodiscard]] std::string PopAttributeOrThrow(char const* name);
        template<class T>
        [[nodiscard]] std::optional<T> PopAttributeMaybe(char const* name);
        template<class T>
        [[nodiscard]] T PopAttributeOrThrow(char const* name);
        template<class T>
        void PopAttributeOrThrow(char const* name, T& out);

        /**
         * Return a string giving a full path to the element.
         * @return
         */
        [[nodiscard]] std::string GetFullPath() const;

        [[nodiscard]] std::string GetPath() const;

        // Return the document this element belongs to
        [[nodiscard]] Document& GetDocument() const;

        // Create a new element as a child of this one, returning it.
        Element AddChild(char const* name);

        // Copy the element and add as child of this one, returning it.
        Element CopyAsChild(Element const&);

        // Delete this element
        void Delete();
        // Delete child of this
        void DeleteChild(Element const& el);

        [[nodiscard]] DyingElement PopChildOrNull(char const* name);
        [[nodiscard]] DyingElement PopChildOrThrow(char const* name);

        void SetAttribute(char const* name, char const* value);
        inline void SetAttribute(char const* name, std::string const& value) {
            SetAttribute(name, value.c_str());
        }
        template <typename T>
        void SetAttribute(char const* name, T&& value);

        operator bool() const;

    private:
        // Until we have deducing this, work around the const/mutable variants.
        template <typename Self, typename F>
        [[nodiscard]] static auto transform_impl(Self&& self, F&& f) {
            using R = std::decay_t<std::invoke_result_t<F, Self>>;
            if (self) {
                return std::optional<R>{std::invoke(std::forward<F>(f), std::forward<Self>(self))};
            } else {
                return std::optional<R>{};
            }
        }
        template <typename Self, typename F>
        [[nodiscard]] static auto and_then_impl(Self&& self, F&& f) {
            using R = std::decay_t<std::invoke_result_t<F, Self>>;
            static_assert(std::disjunction_v<
                    std::is_same<R, std::decay_t<Self>>,
                    util::is_optional<R>
            >);
            if (self) {
                return std::invoke(std::forward<F>(f), std::forward<Self>(self));
            } else {
                return R{};
            }
        }

    public:
        // Returns the result of invocation of f on an Element that does
        // contain a value. Otherwise, returns std::nullopt. f must return
        // Element or a specialisation of std::optional.
        template <typename F>
        [[nodiscard]] auto and_then(F&& f) const {
            return and_then_impl(*this, std::forward<F>(f));
        }
        template <typename F>
        [[nodiscard]] auto and_then(F&& f) {
            return and_then_impl(*this, std::forward<F>(f));
        }

        // Returns a std::optional containing the result of invoking f on
        // an Element that contains a value. Otherwise returns
        // std::nullopt.
        template <typename F>
        [[nodiscard]] auto transform(F&& f) const {
            return transform_impl(*this, std::forward<F>(f));
        }
        template <typename F>
        [[nodiscard]] auto transform(F&& f) {
            return transform_impl(*this, std::forward<F>(f));
        }


    private:

        // For use in templates to hide TinyXML API
        [[nodiscard]] char const* GetAttrOrNull(char const* attr) const;

        // Convert a string to a value
        template <typename T>
        void StringToVal(std::string_view name, char const* s, T& val) const;

        /**
         * Recursive function used by GetPath
         * @param el
         * @param ans
         */
        static void GetPathWorker(const tinyxml2::XMLElement* el, std::ostringstream& ans, bool full);

        // Get the document and use that to create an attribute stream with its factory.
        std::ostringstream MakeAttributeStream() const;

        /**
         * Equality and inequality operators for
         * @param left
         * @param right
         * @return
         */
        friend bool operator==(const Element& left, const Element& right);
        friend bool operator!=(const Element& left, const Element& right);
        friend class AttributeNameIterator;
    };

    // Element which deletes itself on destruction
    // Not copyable obvs
    class DyingElement {
        Element elem;
    public:
        explicit DyingElement(Element e);
        ~DyingElement();
        DyingElement(DyingElement const&) = delete;
        DyingElement& operator=(DyingElement const&) = delete;
        DyingElement(DyingElement&&) = default;
        DyingElement& operator=(DyingElement&&) = default;

        Element& operator*();
        Element* operator->();
        //operator Element&();
        operator bool() const;
    };

    struct AttributeNameIteratorSentinel {};

    class AttributeNameIterator
    {
        Element element;
        tinyxml2::XMLAttribute const* attr = nullptr;

    public:
        using difference_type = std::ptrdiff_t;
        using value_type = char const*;
        using reference = char const*&;
        char const* operator*() const;
        //char const** operator->() const;

        AttributeNameIterator() = default;
        AttributeNameIterator(Element const& e);

        AttributeNameIterator& operator++();
        AttributeNameIterator operator++(int);

        friend bool operator==(const AttributeNameIterator& a, const AttributeNameIterator& b);
        inline friend bool operator==(AttributeNameIterator const& it, AttributeNameIteratorSentinel) {
            return it.attr == nullptr;
        }
    };
    struct AttributeNameIterationRange {
        Element parent;
        [[nodiscard]] AttributeNameIterator begin() const;
        [[nodiscard]] AttributeNameIteratorSentinel end() const;
    };

    // Represent Element::Missing()
    struct ChildIteratorSentinel {};

    /**
     * Want to model the ForwardIterator concept
     */
    class ChildIterator
    {
    protected:
        Element parent = Element::Missing();
        Element current = Element::Missing();

    public:
        using difference_type = std::ptrdiff_t;
        using value_type = Element const;
        using reference = Element const&;
        //using iterator_category = std::forward_iterator_tag;

        Element const& operator*() const;
        Element const* operator->() const;

        inline friend bool operator==(ChildIterator const& it, ChildIteratorSentinel) {
            return it.current == Element::Missing();
        }
    };


    // Iterates over all the child elements.
    class UnnamedChildIterator final : public ChildIterator {
    public:
        UnnamedChildIterator() = default;
        UnnamedChildIterator(const Element& elem);
        UnnamedChildIterator(const Element& elem, const Element& pos);
        UnnamedChildIterator& operator++();
        UnnamedChildIterator operator++(int);
        friend bool operator==(const UnnamedChildIterator& a, const UnnamedChildIterator& b);
    };

    // Iterates over child elements with the given name.
    class NamedChildIterator : public ChildIterator {
    protected:
        std::string name;
    public:
        NamedChildIterator() = default;
        NamedChildIterator(const Element& elem, std::string subElemName);
        NamedChildIterator(const Element& elem, std::string subElemName, const Element& pos);
        NamedChildIterator& operator++();
        NamedChildIterator operator++(int);
        friend bool operator==(const NamedChildIterator& a, const NamedChildIterator& b);
    };

    class NamedChildPopIterator : public NamedChildIterator {
    public:
        using NamedChildIterator::NamedChildIterator;
        NamedChildPopIterator& operator++();
        // Note that after increment the old iterator is garbage as elem is deleted.
        void operator++(int);
    };

    /**
     * Equality comparable
     * @param
     * @return
     */
    //bool operator==(const ChildIterator& a, const ChildIterator& b);

    /**
     * Inequality
     * @param
     * @return
     */
    //bool operator!=(const ChildIterator& a, const ChildIterator& b);

    struct NamedIterationRange {
        Element parent;
        std::string name;

        [[nodiscard]] NamedChildIterator begin() const;
        [[nodiscard]] ChildIteratorSentinel end() const;
    };
    struct UnnamedIterationRange {
        Element parent;

        [[nodiscard]] UnnamedChildIterator begin() const;
        [[nodiscard]] ChildIteratorSentinel end() const;
    };

    struct NamedPopIterationRange {
        Element parent;
        std::string name;

        [[nodiscard]] NamedChildPopIterator begin() const;
        [[nodiscard]] ChildIteratorSentinel end() const;
    };
    /** an abstraction for an XML document
     *
     * this class localises the dependency on an external XML library
     */
    class Document
    {
    public:
        // Create an empty document
        Document();
        /**
         * constructor
         *
         * @param $path
         *   the path to the XML file to be read by this object
         */
        Document(const std::filesystem::path& path);
        template <std::invocable F>
        Document(F&& f) : Document() {
            attr_stream_factory = std::forward<F>(f);
        }
        template <std::invocable F>
        Document(const std::filesystem::path& path, F&& f) : Document(path) {
            attr_stream_factory = std::forward<F>(f);
        }

        /** destructor - needed to avoid polluting all users of this with TinyXML */
        ~Document();

        // User declared destructor => move not implicitly defined
        // Also put these in implementation to avoid tinyxml leaks
        Document(Document&&);
        Document& operator=(Document&&);
        // Rule of 5 recommends explicit deletion of copy, even tho not strictly required
        Document(Document const &) = delete;
        Document& operator=(Document const&) = delete;

        // Deep copy the whole thing
        Document DeepCopy() const;

        Element GetRoot() const;

        // Load some XML from a file
        void LoadFile(const std::filesystem::path& path);

        // Load some XML from a string
        void LoadString(char const* data);
        void LoadString(std::string const& data);

        void SaveFile(const std::filesystem::path& path) const;

        // Create a new element as a child of this one, returning it.
        Element AddChild(char const* name);

        // Get the vector of errors recorded in the document
        inline auto& GetErrors() const {
            return error_list;
        }

        inline void ClearErrors() {
            error_list.clear();
        }

        inline void AddError(std::string descr) {
            error_list.emplace_back(std::move(descr));
        }

    private:
        friend Element;
        friend class DyingElement;

        std::unique_ptr<tinyxml2::XMLDocument> xmlDoc;
        // Note that TinyXML-2 doesn't keep track of
        // the doc filename, so we have to do that.
        std::filesystem::path filename;
        std::function<std::ostringstream()> attr_stream_factory;
        std::vector<std::string> error_list;
    };

    /**
     * Base class for XML errors. Should not be instantiated.
     */
    class XmlError : public hemelb::Exception
    {
    protected:
        //
        /**
         * Should not be instantiated, only its subclasses.
         *
         * Requires the Element where the error occured.
         *
         * @param element
         */
        XmlError();

    public:
        /**
         * Supply an human readable error message.
         * @return
         */
        template<typename T>
        XmlError& operator<<(const T& t)
        {
            return static_cast<XmlError&>(Exception::operator<<(t));
        }

    };

    // Indicate that some XML failed to parse
    class ParseError : public XmlError
    {
    public:
        ParseError(const tinyxml2::XMLDocument*);
    };

    /**
     * Indicate that an element does not have a requested attribute
     */
    class AttributeError : public XmlError
    {
    public:
        AttributeError(const Element& n, std::string_view attr_);

        template<typename T>
        AttributeError& operator<<(const T& t)
        {
            return static_cast<AttributeError&>(XmlError::operator<<(t));
        }
    };

    class DeserialisationError : public XmlError
    {
    public:
        DeserialisationError(const Element& el, std::string_view attrName, std::string_view attrVal);

        template<typename T>
        DeserialisationError& operator<<(const T& t)
        {
            return static_cast<DeserialisationError&> (XmlError::operator<<(t));
        }
    };

    /**
     * Indicate that a requested element does not exist. Should not be used
     * directly; use ChildError, ParentError, or SiblingError.
     */
    class ElementError : public XmlError
    {
    protected:
        ElementError(const Element& el);
    public:
        template<typename T>
        ElementError& operator<<(const T& t)
        {
            return static_cast<ElementError&>(XmlError::operator<<(t));
        }

    protected:
        std::string elemPath;
    };
    /**
     * Indicate that an element lacks the requested child
     */
    class ChildError : public ElementError
    {
    public:
        ChildError(const Element& elem, std::string_view subElemName);
        template<typename T>
        ChildError& operator<<(const T& t)
        {
            return static_cast<ChildError&>(ElementError::operator<<(t));
        }

    };

    /**
     * Indicate that an element lacks a parent.
     */
    class ParentError : public ElementError
    {
    public:
        ParentError(const Element& elem);
        template<typename T>
        ParentError& operator<<(const T& t)
        {
            return static_cast<ParentError&>(ElementError::operator<<(t));
        }

    };

    /**
     * Indicate that an element lacks the requested sibling.
     */
    class SiblingError : public ElementError
    {
    public:
        SiblingError(const Element& elem, std::string_view subElemName);
        template<typename T>
        SiblingError& operator<<(const T& t)
        {
            return static_cast<SiblingError&>(ElementError::operator<<(t));
        }

    };

    template<class T>
    void Element::StringToVal(std::string_view name, char const* s, T& out) const
    {
        const auto N = std::strlen(s);
        if (N == 0)
            throw DeserialisationError(*this, name, s) << "Zero length string for value";

        // So, basically parsing of unsigned values varies across platforms
        // such that on OS X 10.8
        //   unsigned ans
        //   stream = istringstream("-1");
        //   stream >> ans;
        // will NOT set ans but will set stream's failbit.
        //
        // On various Linux boxes, this will set ans == 2**32 - 1 and NOT set
        // the fail bit.
        //
        // Here, for types numeric_limits knows are unsigned we explicitly
        // check for a "-" and throw a suitable error.
        if constexpr (std::numeric_limits<T>::is_specialized && !std::numeric_limits<T>::is_signed)
        {
            if (s[0] == '-')
            {
                throw DeserialisationError(*this, name, s)
                        << " attempt to convert negative number to unsigned type";
            }
        }
        std::istringstream attrStream(s, std::ios_base::in);

        // Don't skip whitespace as that could indicate a malformed value
        // Read bools as true/false
        attrStream >> std::noskipws >> std::boolalpha;

        attrStream >> out;
        if (attrStream.fail())
        {
            throw DeserialisationError(*this, name, s) << " error in extraction operator";
        }
        bool eof = attrStream.eof();
        if (!eof && attrStream.tellg() != std::streampos(N))
        {
            throw DeserialisationError(*this, name, s) << " not all characters consumed";
        }
    }

    template<class T>
    std::optional<T> Element::GetAttributeMaybe(char const* name) const
    {
        const char* attrString = GetAttrOrNull(name);
        if (attrString == nullptr) {
            return std::nullopt;
        } else {
            // Requires T is  default constructible
            std::optional<T> ans{std::in_place};
            StringToVal(name, attrString, *ans);
            return ans;
        }
    }

    template<class T>
    void Element::GetAttributeOrThrow(char const* name, T& out) const
    {
        const char* attrString = GetAttrOrNull(name);
        if (attrString == nullptr)
            throw AttributeError(*this, name);
        StringToVal(name, attrString, out);
    }

    template<class T>
    T Element::GetAttributeOrThrow(char const* name) const {
        T out;
        GetAttributeOrThrow(name, out);
        return out;
    }

    template <typename T>
    std::optional<T> Element::PopAttributeMaybe(const char *name) {
        auto ans = GetAttributeMaybe<T>(name);
        if (ans)
            DelAttribute(name);
        return ans;
    }
    template <typename T>
    T Element::PopAttributeOrThrow(const char *name) {
        auto ans = GetAttributeOrThrow<T>(name);
        DelAttribute(name);
        return ans;
    }

    template<class T>
    void Element::PopAttributeOrThrow(const char *name, T &out) {
        GetAttributeOrThrow(name, out);
        DelAttribute(name);
    }

    template <typename T>
    void Element::SetAttribute(char const* name, T &&value) {
        auto attrStream = MakeAttributeStream();
        attrStream << value;
        SetAttribute(name, attrStream.str().c_str());
    }
}
#endif  // HEMELB_IO_XML_H
